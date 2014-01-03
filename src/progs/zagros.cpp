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

struct DE {
  std::string type;
  size_t position;
  std::string base;
  std::string read;
};

static size_t
convertString(const string& s) {
  std::istringstream buffer(s);
  size_t value;
  buffer >> value;
  return value;
}


static void
sift_single_chrom(vector<GenomicRegion> &other_regions,
                  vector<GenomicRegion> &regions,
                  vector<GenomicRegion> &good_regions) {

  typedef vector<GenomicRegion>::iterator region_itr;

  for (size_t i = 0; i < regions.size(); ++i) {
    region_itr closest(find_closest(other_regions, regions[i]));
    if (closest->overlaps(regions[i])) {
      good_regions.push_back(regions[i]);
      good_regions.back().set_name(closest->get_name());
    }
  }
}

static void
sift(vector<GenomicRegion> &other_regions,
     vector<GenomicRegion> &regions) {

  vector<vector<GenomicRegion> > other_regions_by_chrom;
  separate_chromosomes(other_regions, other_regions_by_chrom);

  unordered_map<string, size_t> chrom_lookup;
  for (size_t i = 0; i < other_regions_by_chrom.size(); ++i)
    chrom_lookup[other_regions_by_chrom[i].front().get_chrom()] = i;
  const vector<GenomicRegion> dummy;

  vector<vector<GenomicRegion> > regions_by_chrom;
  separate_chromosomes(regions, regions_by_chrom);
  regions.clear();

  vector<GenomicRegion> good_regions;
  for (size_t i = 0; i < regions_by_chrom.size(); ++i) {
    const unordered_map<string, size_t>::const_iterator j = 
      chrom_lookup.find(regions_by_chrom[i].front().get_chrom());
    if (j != chrom_lookup.end()) {
      sift_single_chrom(other_regions_by_chrom[j->second], 
			regions_by_chrom[i], good_regions);
    }

  }
  regions.swap(good_regions);

  unordered_map<string, size_t> name_lookup;
  for (size_t i = 0; i < other_regions.size(); ++i)
    name_lookup[other_regions[i].get_name()] = i;

  string prev_name = "";
  for (size_t i = 0; i < regions.size(); ++i)
    if (regions[i].get_name() != prev_name) {
      const size_t idx(name_lookup[regions[i].get_name()]);
      other_regions[idx].set_strand(regions[i].get_strand());
      prev_name = regions[i].get_name();
    }
}

static void
parse_diagnostic_events(const string &de_string,
                        vector<DE> &des) {
  vector<string> parts(smithlab::split_whitespace_quoted(de_string));
  if (parts.size() > 0)
    for (size_t i = 0; i < parts.size(); ++i) {
      DE de;
      size_t index = parts[i].find(">");
      if (index != string::npos) {
        de.type = "mutation";
        de.position = convertString(parts[i].substr(0, index - 1));
        de.base = parts[i].substr(index - 1, 1);
        de.read = parts[i].substr(index + 1, 1);
      } 
      else {
        index = parts[i].find("+");
        if (index != string::npos) {
          de.type = "deletion";
          de.position = convertString(parts[i].substr(0, index));
          de.base = parts[i].substr(index + 1, parts[i].length() - index - 1);
        } 
	else {
          index = parts[i].find("-");
          if (index != string::npos) {
            de.type = "insertion";
            de.position = convertString(parts[i].substr(0, index));
            de.base = parts[i].substr(index + 1, parts[i].length() - index - 1);
          }
        }
      }
      des.push_back(de);
    }
}

static void
add_diagnostic_events_iCLIP(const MappedRead &region,
                            vector<GenomicRegion> &diagnostic_events) {

  vector<DE> des;
  const string name("X");
  parse_diagnostic_events(region.scr, des);
  bool contains_de = true;
  if (des.size() > 0)
    for (size_t i = 0; i < des.size(); ++i)
      if ((region.r.pos_strand() && des[i].type == "insertion"
	   && des[i].base == "T")
          || (!region.r.pos_strand() && des[i].type == "insertion"
              && des[i].base == "A"))
        contains_de = false;
  if (contains_de) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events.back().get_name() == name)
        score = diagnostic_events.back().get_score() + 1;
    if (region.r.pos_strand()) {
      GenomicRegion gr(region.r.get_chrom(), 
		       region.r.get_start(), 
		       region.r.get_start() + 1,
		       name, score, region.r.get_strand());
      diagnostic_events.push_back(gr);
    } 
    else {
      GenomicRegion gr(region.r.get_chrom(), region.r.get_end(), 
		       region.r.get_end() + 1,
		       name, score, region.r.get_strand());
      diagnostic_events.push_back(gr);
    }
  }
}

static void
add_diagnostic_events_hCLIP(const MappedRead &region,
                            vector<GenomicRegion> &diagnostic_events) {
  
  vector<DE> des;
  const string name("X");
  parse_diagnostic_events(region.scr, des);
  if (des.size() == 1) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events.back().get_name() == name)
        score = diagnostic_events.back().get_score() + 1;
    for (size_t j = 0; j < des.size(); ++j) {
      if (region.r.pos_strand() && 
	  des[j].type == "insertion" && des[j].base == "T") {
        GenomicRegion gr(region.r.get_chrom(), 
			 region.r.get_start() + des[j].position - 1,
			 region.r.get_start() + des[j].position, name, score,
			 region.r.get_strand());
        diagnostic_events.push_back(gr);
      } 
      else if (!region.r.pos_strand() && 
	       des[j].type == "insertion" && des[j].base == "A") {
        GenomicRegion gr(region.r.get_chrom(), 
			 region.r.get_start() + des[j].position - 2,
			 region.r.get_start() + des[j].position - 1, name, score,
			 region.r.get_strand());
        diagnostic_events.push_back(gr);
      }
    }
  }
}

static void
add_diagnostic_events_pCLIP(const MappedRead &region,
                            vector<GenomicRegion> &diagnostic_events) {
  vector<DE> des;
  parse_diagnostic_events(region.scr, des);
  const string name("X");
  if (des.size() == 1) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events.back().get_name() == name)
        score = diagnostic_events.back().get_score() + 1;
    for (size_t j = 0; j < des.size(); ++j) {
      if (region.r.pos_strand() && des[j].type == "mutation"
          && des[j].base == "T" && des[j].read == "C") {
        GenomicRegion gr(region.r.get_chrom(), 
			 region.r.get_start() + des[j].position - 1,
			 region.r.get_start() + des[j].position, name, score,
			 region.r.get_strand());
        diagnostic_events.push_back(gr);
      } else if (!region.r.pos_strand() && des[j].type == "mutation"
		 && des[j].base == "A" && des[j].read == "G") {
        GenomicRegion gr(region.r.get_chrom(), 
			 region.r.get_start() + des[j].position - 2,
			 region.r.get_start() + des[j].position - 1, name, score,
			 region.r.get_strand());
        diagnostic_events.push_back(gr);
      }
    }
  }
}

static void
add_diagnostic_events(const MappedRead &region,
                      const string experiment,
                      vector<GenomicRegion> &diagnostic_events) {
  if (experiment == "iCLIP")
    add_diagnostic_events_iCLIP(region, diagnostic_events);
  else if (experiment == "hCLIP")
    add_diagnostic_events_hCLIP(region, diagnostic_events);
  else if (experiment == "pCLIP")
    add_diagnostic_events_pCLIP(region, diagnostic_events);
  else
    throw SMITHLABException("The technology was not recognized!");
}

static void
add_diagnostic_events(const vector<MappedRead> &regions,
                      const string experiment,
                      vector<GenomicRegion> &diagnostic_events) {

  for (size_t i = 0; i < regions.size(); ++i) {
    add_diagnostic_events(regions[i], experiment, diagnostic_events);
  }
}


static void
separate_mapped_reads_chromosomes(const vector<MappedRead> &regions,
                                  vector<vector<MappedRead> > &separated_by_chrom) {
  typedef unordered_map<string, vector<MappedRead> > Separator;
  Separator separator;
  for (vector<MappedRead>::const_iterator i = regions.begin();
       i != regions.end(); ++i) {
    const string the_chrom(i->r.get_chrom());
    if (separator.find(the_chrom) == separator.end())
      separator[the_chrom] = vector<MappedRead>();
    separator[the_chrom].push_back(*i);
  }
  separated_by_chrom.clear();
  for (Separator::iterator i = separator.begin(); i != separator.end(); ++i)
    separated_by_chrom.push_back(i->second);
}

static void
make_inputs(const vector<MappedRead> &mapped_reads,
            const string experiment,
            vector<GenomicRegion> &diagnostic_events) {
  vector<vector<MappedRead> > separated_by_chrom;
  separate_mapped_reads_chromosomes(mapped_reads, separated_by_chrom);
  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    add_diagnostic_events(separated_by_chrom[i], experiment, diagnostic_events);
    separated_by_chrom[i].clear();
  }
}

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
/////////////
/////////////  INPUT DATA TEST AND REFINEMENT
/////////////

static bool
check_overlapping_chrom(const vector<GenomicRegion> &targets_chrom) {

  vector<pair<size_t, bool> > boundaries;
  for (size_t i = 0; i < targets_chrom.size(); ++i) {
    boundaries.push_back(std::make_pair(targets_chrom[i].get_start(), false));
    boundaries.push_back(std::make_pair(targets_chrom[i].get_end(), true));
  }
  sort(boundaries.begin(), boundaries.end());

  size_t count = 0;
  for (size_t i = 0; i < boundaries.size(); ++i)
    if (boundaries[i].second)
      --count;
    else {
      ++count;
      if (count > 1) {
        cerr << targets_chrom[i].tostring() << endl;
        return true;
      }
    }
  return false;
}

static bool
check_overlapping(const vector<GenomicRegion> &targets) {

  vector<vector<GenomicRegion> > targets_chroms;
  separate_chromosomes(targets, targets_chroms);
  bool ret_val = false;
  for (size_t i = 0; i < targets_chroms.size() && !ret_val; ++i)
    ret_val |= check_overlapping_chrom(targets_chroms[i]);
  return ret_val;
}

static char
sample_nuc(const Runif &rng,
           vector<double> &probs) {
  const double d = rng.runif(0.0, 1.0);
  if (d < probs[0])
    return 'A';
  if (d < probs[1])
    return 'C';
  if (d < probs[2])
    return 'G';
  return 'T';
}

static void
replace_Ns(vector<string> &sequences) {
  const Runif rng(std::numeric_limits<int>::max());
  vector<double> probs(vector<double>(smithlab::alphabet_size, 
				      1.0/smithlab::alphabet_size));
  for (size_t i = 0; i < sequences.size(); ++i)
    std::replace(sequences[i].begin(), sequences[i].end(), 'N', 
		 sample_nuc(rng, probs));
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
    if (zoops_i[n] >= 0.8 )
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
      if (zoops_i[n] >= 0.8)
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
    size_t motif_width = 6;
    size_t n_motifs = 1;
    string outfile;
    string chrom_dir = "";
    string structure_file;
    string reads_file;
    string mapper = "rmap";
    string experiment = "iCLIP";
    size_t max_de = std::numeric_limits<size_t>::max();
    double level = std::numeric_limits<double>::max();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<target_regions/sequences>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("width", 'w', "width of motifs to find (4 <= w <= 12; default: 6)", 
		      false, motif_width);
    opt_parse.add_opt("number", 'n', "number of motifs to output (default: 1)", 
		      false, n_motifs);
    opt_parse.add_opt("chrom", 'c', "directory with chrom files (FASTA format)", 
		      false, chrom_dir);
    opt_parse.add_opt("structure", 't', "structure information file", 
		      false, structure_file);
    opt_parse.add_opt("diagnostic_events", 'd', "mapped reads file", 
		      false, reads_file);
    opt_parse.add_opt("level", 'l', "level of influence by diagnostic events (default: maximum)",
          false, level);
    opt_parse.add_opt("mapper", 'm', "Mapper (novoaling, bowtie, rmap; default: rmap)", 
		      false, mapper);
    opt_parse.add_opt("experiment", 'e', 
		      "Type of experiment (hCLIP, pCLIP, iCLIP; default: iCLIP)", 
		      false, experiment);
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
    if (motif_width < 4 || motif_width > 12) {
      cerr << "motif width should be between 4 and 12" << endl;
      return EXIT_SUCCESS;
    }
    if (chrom_dir == "") {
      cerr << "Please specify the directory containing the "
           << "chromosome files" << endl;
      return EXIT_SUCCESS;
    } 
    if (!isdir(chrom_dir.c_str())) {
      throw SMITHLABException(chrom_dir + " is not a valid directory");
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string targets_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    const size_t padding = 0;

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // Data structures and input preparation for sequence
    if (VERBOSE)
      cerr << "LOADING SEQUENCES" << endl;
    vector<string> seqs, names;
    vector<GenomicRegion> targets;
    load_sequences(chrom_dir, padding, targets_file, names, seqs, targets);
    replace_Ns(seqs);

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
    if (!reads_file.empty()) {
      if (VERBOSE)
        cerr << "LOADING MAPPING INFORMATION" << endl;
      vector<MappedRead> mapped_reads;
      load_mapped_reads(reads_file, mapper, mapped_reads);
      if (VERBOSE)
        cerr << "PROCESSING MAPPED READS" << endl;
      if (mapped_reads.size() > 0)
        make_inputs(mapped_reads, experiment, de_regions);
      if (de_regions.size() > 0) {
        sort(de_regions.begin(), de_regions.end());
        if (check_overlapping(targets))
          throw SMITHLABException("Target regions are overlapping!");
        else
          sift(targets, de_regions);
      }
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

//      string de_outfile = outfile.substr(0,outfile.length()-4) + "_de.bed";
//      std::ostream* out3 =
//         (!de_outfile.empty()) ? new std::ofstream(de_outfile.c_str()) : &std::cout;
//      copy(de_regions_sampled.begin(), de_regions_sampled.end(),
//           std::ostream_iterator<GenomicRegion>(*out3, "\n"));
//      if (out3 != &std::cout)
//        delete out3;
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

      model.expectation_maximization(seqs, diagnostic_events,
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
