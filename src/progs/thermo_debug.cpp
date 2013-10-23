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
#include <numeric>

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
  for (size_t i = 0; i < sequences.size(); ++i)
    std::replace(
        sequences[i].begin(), sequences[i].end(), 'N', sample_nuc(rng, probs));
}


int main(int argc, const char **argv) {

  int random_number_seed = numeric_limits<int>::max();
  const Runif rng(random_number_seed);

  try {
    size_t str_level = std::numeric_limits<size_t>::max();

    bool VERBOSE = false;
    string outfile;
    string chrom_dir;

    const size_t structure_computation_window = 20;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(
        strip_path(argv[0]), "", "<sequences>");
    opt_parse.add_opt(
        "output", 'o', "output file name (default: stdout)", false, outfile);
    opt_parse.add_opt("str_level", 'a', "level of structure information",
          false, str_level);
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
    const string targets_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "LOADING SEQUENCES" << endl;
    //Vectors to store primary information from the data
    vector<string> seqs, names;
    vector<vector<double> > secondary_structure;

    read_fasta_file(targets_file, names, seqs);
    replace_Ns(seqs);
    RNAUtils::get_base_pair_probability_vector(VERBOSE, seqs, secondary_structure);

    double level = (100.0 - str_level)/100.0;
    vector<size_t> indices;
    for (size_t i = 0; i < std::floor(level*seqs.size()); ++i) {
      size_t r;
      bool exists = true;
      while (exists) {
        r = rng.runif((size_t)0, secondary_structure.size());
        exists = false;
        for (size_t indx = 0; indx < indices.size(); ++indx)
          if (indices[indx] == r) {
            exists = true;
            break;
          }
        if (!exists)
          indices.push_back(r);
      }
      for (size_t j = 0; j < secondary_structure[r].size(); ++j)
        secondary_structure[r][j] = 0.5;
    }

    save_structure_file(
        secondary_structure, outfile, structure_computation_window);

  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
