/*    
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Emad Bahrami-Samani, Philip J. Uren and Andrew D. Smith
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

#include <cassert>
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

using std::stringstream;
using std::ifstream;
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::pair;
using std::numeric_limits;

/***
 * \brief load the diagnostic events from the given filename. The format of
 *          the file should be one sequence per line, each line is a comma
 *          separated list, where each element is the number of diagnostic
 *          events observed at that location in the given sequence.
 * \param fn            the filename to read the events from
 * \param diagEvents    the events are added to this vector. Any existing data
 *                      is cleared from the vector. diagEvents[i][j] is the
 *                      location of the jth diagnostic event in sequence i.
 *                      locations are relative to the start of the sequence.
 * \return the count of diagnostic events that were found
 * \todo this should be moved into IO.cpp
 */
size_t
loadDiagnosticEvents(const string &fn, vector<vector<size_t> > &diagEvents) {
  size_t total = 0;
  static const size_t buffer_size = 100000; // TODO magic number
  ifstream in(fn.c_str());
  if (!in) throw SMITHLABException("failed to open input file " + fn);

  diagEvents.clear();
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw SMITHLABException("Line too long in file: " + fn);
    vector<string> parts(smithlab::split(string(buffer), ",", false));
    if (parts.size() == 0) continue;
    diagEvents.push_back(vector<size_t>());
    for (size_t j = 0; j < parts.size(); ++j) {
      size_t countAtJ = static_cast<size_t>(atoi(parts[j].c_str()));
      for (size_t a = 0; a < countAtJ; ++a) {
        diagEvents.back().push_back(j);
        total += 1;
      }
    }
  }
  return total;
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

/**
 * \brief TODO
 * \param prob      TODO
 * \param seq_len   TODO
 * \return TODO
 * this is just Poisson probability for 0 observations
 */
static double
prob_no_occurrence(const double prob,
                   const size_t seq_len) {
  return std::exp(-static_cast<int>(seq_len) * prob);
}

/**
 * \brief TODO
 * \param sequences TODO
 * \param base_comp TODO
 */
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

/**
 * \brief TODO
 * \param kmer          TODO
 * \param base_comp     TODO
 * \return TODO
 */
static double
compute_kmer_prob(const string &kmer,
                  const vector<double> &base_comp) {
  double prob = 1.0;
  for (size_t i = 0; i < kmer.length(); ++i)
    prob *= base_comp[base2int(kmer[i])];
  return prob;
}

/**
 * \brief TODO
 * \param kmer         TODO
 * \param base_comp    TODO
 * \param lengths      TODO
 * \return TODO
 */
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

// PJU: following aren't used currently; should they be removed?
/*
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
}*/

// TODO should be in RNA utils
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


/****
 * \TODO this needs to be replaced by use of RNG from smithlab
 */
double
randomDouble(const double lower, const double upper) {
  if (lower >= upper) {
    stringstream ss;
    ss << "Failed to generate random double: "
       << "lower boundary was greater than upper.";
    throw SMITHLABException(ss.str());
  }
  double frac = (double) rand() / RAND_MAX;
  return lower + (frac * (upper - lower));
}


/***
 * \TODO more stuff copied from simulate that needs to go into some library
 */
namespace RNAUT {
  inline size_t
  base2int(char c) {
    switch(c) {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    case 'U' : return 3;
    case 'a' : return 0;
    case 'c' : return 1;
    case 'g' : return 2;
    case 't' : return 3;
    case 'u' : return 3;
    default  : return 4;
    }
  }
}

/****
 * \summary: generate a sequence from a nucleotide distribution.
 * \param dist: the distribution to use.
 * \param length: the length of the sequence to build
 * \throw SMITHLABException: if the dist. vector has the wrong dimensions
 * \TODO this is duplicated in simulate program; it should be pushed into some
 *       common file somewhere. RNA_UTILS probably.
 */
string
genSeqFromNucDist(const vector<double> &dist, const size_t length) {
  if (dist.size() != 4) {   // Todo fix this magic.
    stringstream ss;
    ss << "Failed to generate sequence from nucleotide distribution, "
       << "distribution vector was malformed: found "
       << dist.size() << " entries; expected " << 4; // Todo fix this magic.
    throw SMITHLABException(ss.str());
  }

  string res = "";
  for (size_t i = 0; i < length; ++i) {
    double r = randomDouble(0, 1);
    if (r < dist[RNAUT::base2int('A')]) res += 'A';
    else if (r < dist[RNAUT::base2int('A')] +\
                 dist[RNAUT::base2int('C')]) res += 'C';
    else if (r < dist[RNAUT::base2int('A')] +\
                 dist[RNAUT::base2int('C')] +\
                 dist[RNAUT::base2int('G')]) res += 'G';
    else res += 'T';
  }

  return res;
}


/***
 * \summary given a set of sequences and indicators for motif occurrences,
 *          mask out the most likely occurrences of the motif in the sequences
 * \param seqs          TODO
 * \param indicators    TODO
 * \param zoops         TODO
 * \throw SMITHLABException if the dimensions of seqs doesn't match indicators
 *        or zoops_i
 */
static void
maskOccurrences(vector<string> &seqs, const vector<vector<double> > &indicators,
                const vector<double> &zoops, const size_t motifLen) {
  // todo too much magic here, also this is duplicating stuff in simulate..
  const static double DEFAULT_BACKGROUND[4] = {0.3,  0.2,  0.2,  0.3};
  const static vector<double> DEFAULT_BACKGROUND_VEC
    (DEFAULT_BACKGROUND,
     DEFAULT_BACKGROUND + sizeof(DEFAULT_BACKGROUND) / sizeof(DEFAULT_BACKGROUND[0]));

  if (seqs.size() != indicators.size()) {
    stringstream ss;
    ss << "failed to mask motif occurrences, number of indicator vectors ("
       << indicators.size() << ") didn't match the number of sequences ("
       << seqs.size() << ")";
    throw SMITHLABException(ss.str());
  }
  if (seqs.size() != zoops.size()) {
    stringstream ss;
    ss << "failed to mask motif occurrences, number of zoops indicators ("
       << zoops.size() << ") didn't match the number of sequences ("
       << seqs.size() << ")";
    throw SMITHLABException(ss.str());
  }


  for (size_t i = 0; i < seqs.size(); ++i) {
    if (zoops[i] >= 0.8) {      // TODO Abracadabra!
      double max_X = -1;
      int max_indx = -1;
      for (size_t j = 0; j < indicators[i].size(); j++) {
        if (indicators[i][j] > max_X) {
          max_X = indicators[i][j];
          max_indx = j;
        }
      }
      string junk = genSeqFromNucDist(DEFAULT_BACKGROUND_VEC, motifLen);
      seqs[i].replace(max_indx, motifLen, junk);
    }
  }
}

static string
format_motif_header(const string &name) {
  static const string the_rest("XX\nTY\tMotif\nXX\nP0\tA\tC\tG\tT");
  std::ostringstream oss;
  oss << "AC\t" << name << '\n' << the_rest;
  return oss.str();
}

/***
 * \brief TODO
 * \param model         TODO
 * \param motif_name    TODO
 * \param targets       TODO
 * \param sequences     TODO
 * \param indicators    TODO
 * \param zoops_i       TODO
 * \return  TODO
 * \throw   TODO
 */
static string
format_motif(const Model &model,
             const string &motif_name,
             const vector<GenomicRegion> &targets,
             const vector<string> &sequences,
             const vector<vector<double> > &indicators,
             const vector<double> &zoops_i) {
  if (sequences.size() != indicators.size()) {
    stringstream ss;
    ss << "failed to format motif for output. Number of site indicator vectors "
       << "(" << indicators.size() << ") does not equal number of sequences ("
       << sequences.size() << ")";
    throw SMITHLABException(ss.str());
  }
  if (zoops_i.size() != sequences.size()) {
    stringstream ss;
    ss << "failed to format motif for output. Number of sequence indicators "
       << "(" << zoops_i.size() << ") " << "does not match number of sequences "
       << sequences.size();
    throw SMITHLABException(ss.str());
  }

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
    // TODO fix below
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
    if (!targets.empty()) {
      if (zoops_i[n] >= 0.8)
        ss << format_site(model, targets[n], sequences[n], site_pos) << "\t"
	       << zoops_i[n] << endl;
    }
  }
  ss << "XX" << endl << "//";
  return ss.str();
}

int main(int argc, const char **argv) {
  try {
    static const double zoops_expansion_factor = 0.75;

    // TODO: the level parameter, specifying how much influence the DEs have,
    // is currently not used
    bool VERBOSE = false;
    size_t motif_width = 6;
    size_t n_motifs = 1;
    string outfile;
    string chrom_dir = "";
    string structure_file;
    string reads_file;
    // not implemented
    //double level = std::numeric_limits<double>::max();
    size_t numStartingPoints = 3;

    /****************** COMMAND LINE OPTIONS ********************/
    // TODO -- specify any missing default values below using constants and
    // fix the ones that are hard-coded in.
    OptionParser opt_parse(strip_path(argv[0]), "", "<target_regions/sequences>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)", 
                      OptionParser::OPTIONAL, outfile);
    opt_parse.add_opt("width", 'w', "width of motifs to find (4 <= w <= 12; "
                      "default: 6)", OptionParser::OPTIONAL, motif_width);
    opt_parse.add_opt("number", 'n', "number of motifs to output (default: 1)", 
                      OptionParser::OPTIONAL, n_motifs);
    opt_parse.add_opt("chrom", 'c', "directory with chrom files (FASTA format)", 
                      OptionParser::OPTIONAL, chrom_dir);
    opt_parse.add_opt("structure", 't', "structure information file", 
                      OptionParser::OPTIONAL, structure_file);
    opt_parse.add_opt("diagnostic_events", 'd', "mapped reads file", 
                      OptionParser::OPTIONAL, reads_file);
    // not currently implemented
    /*opt_parse.add_opt("level", 'l', "level of influence by diagnostic events "
                      "(default: maximum)", OptionParser::OPTIONAL, level);*/
    opt_parse.add_opt("starting-points", 's', "number of starting points to try "
                      "for EM search. Higher values will be slower, but more "
                      "likely to find the global maximum.",
                      OptionParser::OPTIONAL, numStartingPoints);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      OptionParser::OPTIONAL, VERBOSE);
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
    if (leftover_args.size() != 1) {
      cerr << "Zagros requires one input file, found "
           << leftover_args.size() << ": ";
      for (size_t i = 0; i < leftover_args.size(); ++i) {
        cerr << leftover_args[i];
        if (i != (leftover_args.size() - 1)) cerr << ", ";
      }
      cerr << endl << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string targets_file(leftover_args.back());
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
    load_sequences(targets_file, chrom_dir, seqs, names, targets);
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

    // Load the diagnostic events
    vector<vector<size_t> > diagEvents (seqs.size());
    if (!reads_file.empty()) {
      if (VERBOSE)
        cerr << "LOADING DIAGNOSTIC EVENTS... ";
      size_t deCount = loadDiagnosticEvents(reads_file, diagEvents);
      if (diagEvents.size() != seqs.size()) {
        stringstream ss;
        ss << "inconsistent dimensions of sequence and diagnostic events data. "
           << "Found " << seqs.size() << " sequences, and " << diagEvents.size()
           << " diagnostic events vectors";
        throw SMITHLABException(ss.str());
      }
      if (VERBOSE)
        cerr << "DONE (FOUND " << deCount << " EVENTS IN TOTAL)" << endl;
    }

    // find and output each of the motifs that the user asked for.
    for (size_t i = 0; i < n_motifs ; ++i) {
      if (VERBOSE)
        cerr << "FITTING MOTIF PARAMETERS FOR MOTIF " << (i+1)
             << " OF " << n_motifs << endl;

      // pick a set of starting points to try
      vector<kmer_info> top_kmers;
      find_best_kmers(motif_width, numStartingPoints, seqs, top_kmers);
      double bestLogLike = 0;
      bool firstKmer = true;

      vector<vector<double> > indicators;
      vector<double> has_motif;
      Model model;

      for (size_t j = 0; j < numStartingPoints; ++j) {
        Model model_l;
        Model::set_model_by_word(Model::pseudocount, top_kmers[j].kmer, model_l);
        model_l.p = 0.5;
        model_l.gamma = ((seqs.size() - (zoops_expansion_factor*
                       (seqs.size() - top_kmers[j].observed)))/
               static_cast<double>(seqs.size()));
        if (!secondary_structure.empty()) {
          model_l.motif_sec_str = vector<double>(motif_width, 0.5);
          model_l.f_sec_str = 0.5;
        }

        vector<double> has_motif_l(seqs.size(), model_l.gamma);
        vector<vector<double> > indicators_l;

        for (size_t k = 0; k < seqs.size(); ++k) {
          const size_t n_pos = seqs[k].length() - motif_width + 1;
          indicators_l.push_back(vector<double>(n_pos, 1.0 / n_pos));
        }

        if (VERBOSE)
          cerr << "\t" << "TRYING STARTING POINT " << (j+1) << " OF "
          << numStartingPoints << " (" << top_kmers[j].kmer << ") ... ";
        model_l.expectationMax(seqs, diagEvents, secondary_structure,
            indicators_l, has_motif_l);
        double logLike;
        if (secondary_structure.size() == 0) {
          logLike = model_l.calculate_zoops_log_l(seqs, diagEvents,
                                                  indicators_l, has_motif_l);
        } else {
          logLike = model_l.calculate_zoops_log_l(seqs, secondary_structure,
                                                  diagEvents, indicators_l,
                                                  has_motif_l);
        }
        if (VERBOSE)
          cerr << "LOG-LIKELIHOOD: " << logLike << endl;

        if ((firstKmer) || (logLike > bestLogLike)) {
          bestLogLike = logLike;
          model = model_l;
          indicators = indicators_l;
          has_motif = has_motif_l;
          firstKmer = false;
        }
      }

      if (VERBOSE)
        cerr << "\t" << "WRITING MOTIF " << endl;
      out << format_motif(model, "ZAGROS" + toa(i), targets, seqs,
			  indicators, has_motif) << endl;

      if (VERBOSE)
              cerr << "\t" << "MASKING MOTIF OCCURRENCES" << endl;
      maskOccurrences(seqs, indicators, has_motif, motif_width);
    }
  } 
  catch (const SMITHLABException &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  } 
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
