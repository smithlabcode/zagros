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
#include <sstream>
#include <iterator>
#include <fstream>
#include <numeric>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "Model.hpp"
#include "IO.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ofstream;
using std::ostream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::sort;
using std::accumulate;
using std::numeric_limits;

using smithlab::alphabet_size;

int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outdir = "zagros_output";
    string chrom_dir;
    size_t motif_width = 8;
    string experiment = "none";
    string mapper;
    string structure_information_file = "";
    bool use_de_information = false;

    string reads_file;

    const size_t max_iterations = 10;
    const double tolerance = 1e-10;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(
        strip_path(argv[0]), "", "<target regions/sequences>");
    opt_parse.add_opt(
        "output", 'o', "Name of output directory (default: zagros_output)",
        false, outdir);
    opt_parse.add_opt("width", 'w', "motif width", false, motif_width);
    opt_parse.add_opt(
        "chrom", 'c', "directory with chrom files (FASTA format)", false,
        chrom_dir);
//    opt_parse.add_opt(
//        "experiment",
//        'e',
//        "The type of experiment: iCLIP, hCLIP, pCLIP (default: no diagnostic events considered)",
//        false, experiment);
//    opt_parse.add_opt("mapper", 'm', "Mapped reads format: novoalign, "
//        "bowtie, rmap, piranha (default: none)", false, mapper);
//    opt_parse.add_opt("reads", 'r', "Mapped reads file", false, reads_file);
//    opt_parse.add_opt(
//        "structure", 't', "Use the structure information", false,
//        structure_information_file);
//    opt_parse.add_opt(
//        "diagnostic_events", 'd', "Use the diagnostic events information",
//        false, use_de_information);
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

    //Create the output for directory
    if (outdir.find_last_of("/\\") != outdir.length() - 1)
      outdir += "/";
    mkdir(outdir.c_str(), 0750);

    string dirname, base_name, suffix;
    parse_dir_baseanme_suffix(targets_file, dirname, base_name, suffix);
    string base_file = outdir + base_name;

    //Vectors to store primary information from the data
    vector<string> seqs;
    vector<string> names;
    vector<GenomicRegion> targets;
    vector<vector<size_t> > diagnostic_events;
    vector<double> secondary_structure;

    size_t padding = 0;
    if (structure_information_file != "") {
      std::ifstream in(structure_information_file.c_str(), std::ios::binary);
      if (!in) {
        throw SMITHLABException(
            "cannot open input file " + string(structure_information_file));
      } else {
        in.close();
        padding = IO::flanking_regions_size;
      }
    }

    IO::load_sequences(names, seqs, targets, chrom_dir, padding, targets_file);

    vector<double> has_motif(seqs.size(), 1.0);
    vector<vector<double> > indicators;
    for (size_t i = 0; i < seqs.size(); ++i) {
      size_t window_size =
          (targets.empty()) ? seqs[i].length() : targets[i].get_width();
      const size_t n_pos = window_size - motif_width + 1;
      indicators.push_back(vector<double>(n_pos, 1.0 / n_pos));
    }

    Model model(motif_width);
    model.expectation_maximization(
        seqs, diagnostic_events, secondary_structure, indicators, has_motif,
        structure_information_file, use_de_information, base_file,
        max_iterations, tolerance);

    string output_model = base_file + ".mat";
    std::ostream* outf =
        (!output_model.empty()) ? new ofstream(output_model.c_str()) : &cout;
    *outf
        << IO::print_model(
            model, "ME_EM", targets, seqs, indicators, has_motif);
    delete outf;

  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
