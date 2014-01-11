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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "Part_Func.hpp"
#include "RNA_Utils.hpp"
#include "IO.hpp"
#include "RNG.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::numeric_limits;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////
/////////////  REPLACING THE Ns IN SEQUENCING RANDOMLY
/////////////

static char sample_nuc(const Runif &rng, vector<double> &probs) {
  const double d = rng.runif(0.0, 1.0);
  if (d < probs[0])
    return 'A';
  if (d < probs[1])
    return 'C';
  if (d < probs[2])
    return 'G';
  return 'T';
}

static void replace_Ns(vector<string> &sequences) {
  const Runif rng(std::numeric_limits<int>::max());
  vector<double> probs(
      vector<double>(smithlab::alphabet_size, 1.0 / smithlab::alphabet_size));
  for (size_t i = 0; i < sequences.size(); ++i) {
    for (size_t j = 0; j < sequences[i].size(); j++) {
      if (sequences[i][j] == 'N') {
        std::replace(sequences[i].begin() + j, sequences[i].begin() + j + 1,
                     'N', sample_nuc(rng, probs));
      }
    }
  }
}


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    string chrom_dir;

    const size_t sequencePadding = 20;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(
        strip_path(argv[0]), "", "<target_regions/sequences>");
    opt_parse.add_opt(
        "output", 'o', "output file name (default: stdout)", false, outfile);
    opt_parse.add_opt(
        "chrom", 'c', "directory with chrom files (FASTA format)", false,
        chrom_dir);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string targFile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    // where to write otuput
    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    //Vectors to store primary information from the data
    vector<string> seqs, names;
    vector<GenomicRegion> targets;
    vector<vector<double> > secStructure;

    if (VERBOSE) cerr << "LOADING SEQUENCES... ";
    load_sequences(targFile, chrom_dir, seqs, names, targets, sequencePadding);
    replace_Ns(seqs);
    if (VERBOSE) cerr << "DONE" << endl;

    RNAUtils::get_base_pair_probability_vector(VERBOSE, seqs, secStructure);
    save_structure_file(secStructure, out, sequencePadding);

  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
