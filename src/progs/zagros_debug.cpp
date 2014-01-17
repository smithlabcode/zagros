/*    
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Emad Bahrami-Samani and Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <numeric>
#include <queue>
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "RNG.hpp"

#include "Model.hpp"
#include "IO.hpp"
using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::pair;
using std::numeric_limits;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////
/////////////  STUFF FOR DIAGNOSTIC EVENTS
/////////////

static void
load_diagnostic_events(const vector<GenomicRegion> &regions,
                       const vector<GenomicRegion> &de_regions,
                       const size_t max_de,
                       vector<vector<size_t> > &D) {

  unordered_map<string, size_t> name_lookup;
  for (size_t i = 0; i < regions.size(); ++i)
    name_lookup[regions[i].get_name()] = i;

  D.resize(regions.size());
  for (size_t i = 0; i < de_regions.size(); ++i) {
    const size_t idx = name_lookup[de_regions[i].get_name()];
    if (D[idx].size() < max_de) {
      if (regions[idx].pos_strand())
        D[idx].push_back(de_regions[i].get_start() -
             regions[idx].get_start() - 1);
      else D[idx].push_back(regions[idx].get_width() -
           (de_regions[i].get_start() -
            regions[idx].get_start() + 1));
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////
/////////////  STUFF FOR FINDING THE STARTING POINT BELOW HERE
/////////////

struct kmer_info {
  std::string kmer;
  double expected;
  size_t observed;
  kmer_info(const std::string &km, const double ex, const double ob) :
    kmer(km), expected(ex), observed(ob) {}
  double score() const {
    return observed/expected;
  }
  bool operator>(const kmer_info &ki) const {
    return score() > ki.score();
  }
};

// this is just Poisson probability for 0 observations
static double
prob_no_occurrence(const double prob,
                   const size_t seq_len) {
  return std::exp(-static_cast<int>(seq_len) * prob);
}

static void
compute_base_comp(const vector<string> &sequences,
                  vector<double> &base_comp) {
  base_comp.resize(smithlab::alphabet_size, 0.0);
  size_t total = 0;
  for (size_t i = 0; i < sequences.size(); ++i) {
    for (size_t j = 0; j < sequences[i].length(); ++j)
      ++base_comp[base2int(sequences[i][j])];
    total += sequences[i].length();
  }
  std::transform(
     base_comp.begin(), base_comp.end(), base_comp.begin(),
     std::bind2nd(std::divides<double>(), total));
}

static double
compute_kmer_prob(const string &kmer,
                  const vector<double> &base_comp) {
  double prob = 1.0;
  for (size_t i = 0; i < kmer.length(); ++i)
    prob *= base_comp[base2int(kmer[i])];
  return prob;
}

static double
expected_seqs_with_kmer(const string &kmer,
                        const vector<double> &base_comp,
                        const vector<size_t> &lengths) {
  const double p = compute_kmer_prob(kmer, base_comp);
  double expected = 0.0;
  for (size_t i = 0; i < lengths.size(); ++i)
    expected += (1.0 - prob_no_occurrence(p, lengths[i]));
  return expected;
}

static size_t count_seqs_with_kmer(const string &kmer,
                                   const vector<string> &sequences) {
  size_t count = 0;
  for (size_t i = 0; i < sequences.size(); ++i) {
    bool has_kmer = false;
    const size_t lim = sequences[i].length() - kmer.length() + 1;
    for (size_t j = 0; j < lim && !has_kmer; ++j)
      has_kmer = !sequences[i].compare(j, kmer.length(), kmer);
    count += has_kmer;
  }
  return count;
}

static void
find_best_kmers(const size_t k_value,
                const size_t n_top_kmers,
                const vector<string> &sequences,
                vector<kmer_info> &top_kmers) {

  const size_t n_kmers = (1ul << 2 * k_value);

  vector<double> base_comp;
  compute_base_comp(sequences, base_comp);

  vector<size_t> lengths;
  for (size_t i = 0; i < sequences.size(); ++i)
    lengths.push_back(sequences[i].length());

  std::priority_queue<kmer_info, vector<kmer_info>,
          std::greater<kmer_info> > best_kmers;

  for (size_t i = 0; i < n_kmers; ++i) {
    const string kmer(i2mer(k_value, i));
    const double expected = expected_seqs_with_kmer(kmer, base_comp, lengths);
    const size_t observed = count_seqs_with_kmer(kmer, sequences);
    best_kmers.push(kmer_info(kmer_info(kmer, expected, observed)));
    if (best_kmers.size() > n_top_kmers)
      best_kmers.pop();
  }

  while (!best_kmers.empty()) {
    top_kmers.push_back(best_kmers.top());
    best_kmers.pop();
  }
  reverse(top_kmers.begin(), top_kmers.end());
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////
///////////  CODE FOR FORMATTING MOTIF / MODEL OUTPUT BELOW HERE
///////////

static string
format_site(const Model &model,
            const GenomicRegion &region,
            const string &seq,
            const size_t site_pos) {
  const size_t width = model.size();
  std::ostringstream ss;
  ss << "BS\t" << seq.substr(site_pos, width) << "; "
     << assemble_region_name(region) << "; " << site_pos << "; " << width
     << ";  ;" << (region.pos_strand() ? "p;" : "n;");
  return ss.str();
}

static string
format_motif_header(const string &name) {
  static const string the_rest("XX\nTY\tMotif\nXX\nP0\tA\tC\tG\tT");
  std::ostringstream oss;
  oss << "AC\t" << name << '\n' << the_rest;
  return oss.str();
}

static string
format_motif(const Model &model,
             const string &motif_name,
             const vector<GenomicRegion> &targets,
             const vector<string> &sequences,
             const vector<vector<double> > &indicators,
             const vector<double> &zoops_i) {

  assert(
   sequences.size() == indicators.size()
   && zoops_i.size() == sequences.size());

  std::ostringstream ss;
  ss << format_motif_header(motif_name) << endl;

  vector<vector<double> > tmp_m = model.matrix;
  for (size_t n = 0; n < sequences.size(); n++) {
    double max_X = -1;
    int max_i = -1;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < model.size(); ++j)
      tmp_m[j][base2int(sequences[n][max_i + j])] += zoops_i[n];
  }

  for (size_t j = 0; j < tmp_m.size(); j++) {
    // AS: this is crazy below and will not work for motifs of width 11
    ss << "0" << j + 1;
    for (size_t b = 0; b < smithlab::alphabet_size; ++b)
      ss << '\t' << static_cast<int>(tmp_m[j][b]);
    ss << endl;
  }

  if (model.motif_sec_str.size() > 0) {
    ss << "XX" << endl << "AT\tSEC_STR=";
    for (size_t i = 0; i < model.motif_sec_str.size() - 1; ++i)
      ss << model.motif_sec_str[i] << ",";
    ss << model.motif_sec_str[model.motif_sec_str.size() - 1] << endl;
  }

  ss << "XX" << endl << "AT\tGEO_P=" << model.p << endl;
  ss << "AT\tGEO_DELTA=" << model.delta << endl << "XX" << endl;

  for (size_t n = 0; n < indicators.size(); ++n) {
    double max_X = -1;
    size_t site_pos = 0;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        site_pos = i;
      }
    }
    if (!targets.empty())
      if (zoops_i[n] > Model::zoops_threshold)
        ss << format_site(model, targets[n], sequences[n], site_pos) << "\t"
     << zoops_i[n] << endl;
  }
  ss << "XX" << endl << "//";

  return ss.str();
}

int main(int argc,
         const char **argv) {

  try {

    static const double zoops_expansion_factor = 0.75;

    bool VERBOSE = false;
    size_t motif_width = 8;
    size_t n_motifs = 1;
    string outfile;
    string structure_file;
    size_t max_de = std::numeric_limits<size_t>::max();
    size_t level = std::numeric_limits<size_t>::max();
    size_t de_l = std::numeric_limits<size_t>::max();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<target_regions/sequences>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
          false, outfile);
    opt_parse.add_opt("width", 'w', "width of motifs to find",
          false, motif_width);
    opt_parse.add_opt("number", 'n', "number of motifs to output",
          false, n_motifs);
    opt_parse.add_opt("structure", 't', "structure information file",
          false, structure_file);
    opt_parse.add_opt("de_level", 'a', "level of diagnostic events",
          false, de_l);
    opt_parse.add_opt("level", 'l', "level of influence by diagnostic events",
          false, level);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl << opt_parse.about_message()
     << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string sequence_file(leftover_args.front());
    const string targets_file(leftover_args[1]);
    const string de_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // Data structures and input preparation for sequence
    if (VERBOSE)
      cerr << "LOADING SEQUENCES" << endl;
    vector<string> seqs, names;
    vector<GenomicRegion> targets;
    ReadBEDFile(targets_file, targets);
    read_fasta_file(sequence_file, names, seqs);

    // Data structures and input preparation for secondary structure
    vector<vector<double> > secondary_structure;
    if (!structure_file.empty()) {
      if (VERBOSE)
        cerr << "LOADING STRUCTURE INFORMATION" << endl;
      load_structures(structure_file, secondary_structure);
      if (!seq_and_structure_are_consistent(seqs, secondary_structure))
        throw SMITHLABException("inconsistent dimensions of "
        "sequence and structure data");
    }

    // Data structures and input preparation for diagnostic events
    vector<GenomicRegion> de_regions;
    vector<vector<size_t> > diagnostic_events;
    if (VERBOSE)
      cerr << "LOADING MAPPING INFORMATION" << endl;
    ReadBEDFile(de_file, de_regions);

    int random_number_seed = numeric_limits<int>::max();
    const Runif rng(random_number_seed);

    vector<GenomicRegion> de_regions_sampled;
    if (de_regions.size() > std::floor(level * targets.size()))
      for (size_t i = 0; i < std::floor(level*targets.size()); ++i) {
        const size_t r = rng.runif((size_t)0, de_regions.size());
        de_regions_sampled.push_back(de_regions[r]);
        de_regions.erase(de_regions.begin() + r);
      }
    else
      de_regions_sampled = de_regions;

//      std::ofstream de_out(string(outfile.substr(0,outfile.length()-4) + "_de.bed").c_str());
//      copy(
//          de_regions_sampled.begin(),
//          de_regions_sampled.end(),
//          std::ostream_iterator<GenomicRegion>(de_out, "\n"));
    if (level > 0)
      load_diagnostic_events(targets, de_regions_sampled, max_de, diagnostic_events);

    double de_level = (100.0 - de_l)/100.0;
    vector<size_t> indices;
    for (size_t i = 0; i < std::floor(de_level*seqs.size()); ++i) {
      size_t r;
      bool exists = true;
      while (exists) {
        r = rng.runif((size_t)0, diagnostic_events.size());
        exists = false;
        for (size_t indx = 0; indx < indices.size(); ++indx)
          if (indices[indx] == r) {
            exists = true;
            break;
          }
        if (!exists)
          indices.push_back(r);
      }
      diagnostic_events[r].clear();
    }

    if (VERBOSE)
      cerr << "IDENTIFYING STARTING POINTS" << endl;
    vector<kmer_info> top_kmers;
    find_best_kmers(motif_width, n_motifs, seqs, top_kmers);

    if (VERBOSE)
      cerr << "FITTING MOTIF PARAMETERS" << endl;

    for (size_t i = 0; i < top_kmers.size() ; ++i) {

      Model model;
      Model::set_model_by_word(Model::pseudocount, top_kmers[i].kmer, model);
      model.p = 0.5;
      model.gamma = ((seqs.size() - (zoops_expansion_factor*
             (seqs.size() - top_kmers[i].observed)))/
         static_cast<double>(seqs.size()));
      if (!secondary_structure.empty()) {
        model.motif_sec_str = vector<double>(motif_width, 0.5);
        model.f_sec_str = 0.5;
      }

      vector<double> has_motif(seqs.size(), model.gamma);
      vector<vector<double> > indicators;
      for (size_t j = 0; j < seqs.size(); ++j) {
        const size_t n_pos = seqs[j].length() - motif_width + 1;
        indicators.push_back(vector<double>(n_pos, 1.0 / n_pos));
      }

      model.expectationMax(seqs, diagnostic_events,
             secondary_structure, indicators, has_motif);

      out << format_motif(model, "ZAGROS" + toa(i), targets, seqs,
        indicators, has_motif) << endl;
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
