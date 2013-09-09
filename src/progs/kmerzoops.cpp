/*    kmerzoops: find most enriched k-mers assuming zoops
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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
#include <queue>

#include <OptionParser.hpp>
#include <smithlab_utils.hpp>
#include <smithlab_os.hpp>

using std::string;
using std::vector;
using std::cerr;
using std::endl;



static double
compute_kmer_prob(const string &kmer, const vector<double> &base_comp) {
  double prob = 1.0;
  for (size_t i = 0; i < kmer.length(); ++i)
    prob *= base_comp[base2int(kmer[i])];
  return prob;
}



// this is just Poisson probability for 0 observations
static double
prob_no_occurrence(const double prob, const size_t seq_len) {
  return std::exp(-static_cast<int>(seq_len)*prob);
}



static double
expected_seqs_with_kmer(const string &kmer, const vector<double> &base_comp,
			const vector<size_t> &lengths) {
  const double p = compute_kmer_prob(kmer, base_comp);
  double expected = 0.0;
  for (size_t i = 0; i < lengths.size(); ++i)
    expected += (1.0 - prob_no_occurrence(p, lengths[i]));
  return expected;
}



static size_t
count_seqs_with_kmer(const string &kmer, const vector<string> &sequences) {
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
compute_base_comp(const vector<string> &sequences, vector<double> &base_comp) {
  base_comp.resize(smithlab::alphabet_size, 0.0);
  size_t total = 0;
  for (size_t i = 0; i < sequences.size(); ++i) {
    for (size_t j = 0; j < sequences[i].length(); ++j)
      ++base_comp[base2int(sequences[i][j])];
    total += sequences[i].length();
  }
  std::transform(base_comp.begin(), base_comp.end(), base_comp.begin(), 
  		 std::bind2nd(std::divides<double>(), total));
}



struct kmer_info {
  string kmer;
  double expected;
  size_t observed;
  kmer_info(const string &km, const double ex, const double ob) :
    kmer(km), expected(ex), observed(ob) {}
  double score() const {return observed/expected;}
  bool operator>(const kmer_info &ki) const {return score() > ki.score();}
};



static std::ostream&
operator<<(std::ostream &os, const kmer_info &ki) {
  return os << ki.kmer << '\t' << ki.observed << '\t'
	    << ki.expected << '\t' << ki.score();
}



int 
main(int argc, const char **argv) {
  
  try {

    bool VERBOSE = false;
    string outfile;
    size_t n_top_kmers = 1;
    size_t k_value = 6;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
			   "find most enriched k-mers assuming zoops",
			   "<sequences>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("number", 'n', "number of top k-mers to output",
		      false, n_top_kmers);
    opt_parse.add_opt("kmer", 'k', "width of k-mers", false, k_value);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string sequences_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream in(sequences_file.c_str());
    if (!in)
      throw SMITHLABException("cannot open input file: " + sequences_file);

    vector<string> names, sequences;
    string buffer;
    while (getline(in, buffer)) {
      if (buffer.length() > 0) {
	if (buffer[0] == '>') {
	  names.push_back(buffer.substr(1));
	  sequences.push_back("");
	}
	else if (!sequences.empty())
	  sequences.back() += buffer;
      }
    }

    const size_t n_kmers = (1ul << 2*k_value);
    if (VERBOSE)
      cerr << "N SEQUENCES:\t" << sequences.size() << endl
	   << "K MERS TO CHECK:\t" << n_kmers << endl;

    vector<double> base_comp;
    compute_base_comp(sequences, base_comp);

    if (VERBOSE) {
      cerr << "BASE COMP" << endl;
      for (size_t i = 0; i < base_comp.size(); ++i)
	cerr << int2base(i) << '\t' << base_comp[i] << endl;
    }

    vector<size_t> lengths;
    for (size_t i = 0; i < sequences.size(); ++i)
      lengths.push_back(sequences[i].length());

    std::priority_queue<kmer_info,
			vector<kmer_info>,
			std::greater<kmer_info> > best_kmers;

    if (VERBOSE)
      cerr << "EVALUATING K-MERS" << endl;

    for (size_t i = 0; i < n_kmers; ++i) {
      const string kmer(i2mer(k_value, i));
      const double expected = expected_seqs_with_kmer(kmer, base_comp, lengths);
      const size_t observed = count_seqs_with_kmer(kmer, sequences);
      best_kmers.push(kmer_info(kmer_info(kmer, expected, observed)));
      if (best_kmers.size() > n_top_kmers)
	best_kmers.pop();
    }

    vector<kmer_info> to_reverse;
    while (!best_kmers.empty()) {
      to_reverse.push_back(best_kmers.top());
      best_kmers.pop();
    }
    reverse(to_reverse.begin(), to_reverse.end());
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    copy(to_reverse.begin(), to_reverse.end(), 
	 std::ostream_iterator<kmer_info>(out, "\n"));
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
