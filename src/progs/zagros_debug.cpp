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
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::accumulate;

using smithlab::alphabet_size;

int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    size_t motif_width = 8;

    const size_t max_iterations = 10;
    const double tolerance = 1e-10;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<seqs>");
    opt_parse.add_opt(
        "output", 'o', "Name of output file (default: stdout)", false, outfile);
    opt_parse.add_opt("width", 'w', "motif width", false, motif_width);
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
    if (leftover_args.size() != 4) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string seqs_file(leftover_args.front());
    const string regions_file(leftover_args[1]);
    const string de_file(leftover_args[2]);
    const string motif_no(leftover_args.back());

    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> seqs;
    vector<string> names;
    vector<GenomicRegion> regions;
    vector<GenomicRegion> de_regions;

    read_fasta_file(seqs_file, names, seqs);
    ReadBEDFile(regions_file, regions);
    ReadBEDFile(de_file, de_regions);

    unordered_map<string, size_t> names_table;
    IO::make_sequence_names(names, seqs, regions, names_table);

    vector<vector<size_t> > diagnostic_events;
    IO::load_diagnostic_events(
        de_regions, names_table, regions, diagnostic_events);

    stringstream convert(motif_no);
    size_t mn;
    convert >> mn;
    const string motif_name = "DE_EM_" + motif_no;
    size_t base_index = regions_file.find_last_of("/\\");
    string base_file = regions_file.substr(
        base_index + 1, regions_file.length() - base_index - 5);

    vector<vector<size_t> > tmp_de = diagnostic_events;
    vector<vector<size_t> > tmp_de_c(seqs.size());

    size_t k = 0;
    for (size_t i = 0; i < seqs.size() && k < mn; ++i) {
      for (size_t j = 0; j < diagnostic_events[i].size() && k < mn; j += 1) {
        tmp_de_c[i].push_back(tmp_de[i][j]);
        k++;
      }
    }

    vector<vector<double> > indicators;
    for (size_t ik = 0; ik < seqs.size(); ++ik) {
      const size_t n_pos = seqs[ik].length() - motif_width + 1;
      indicators.push_back(vector<double>(n_pos, 1.0 / n_pos));
    }

    Model model(motif_width, regions);
//    model.find_delta(seqs, regions, diagnostic_events);

    for (size_t i = 0; i<diagnostic_events.size(); ++i) {
      for (size_t j = 0; j<diagnostic_events[i].size(); ++j)
        cout << tmp_de_c[i][j] << "\t";
      cout << endl;
    }

    model.expectation_maximization(
        max_iterations, tolerance, regions, seqs, tmp_de_c, indicators, 1, 0,
        10000, base_file);

//    IO::save_input_files(seqs, regions, de_regions, base_file);
    model.prepare_output(seqs, indicators, tmp_de_c, base_file);

    cout << model.print_model("DE_EM", regions, seqs, indicators);

//      set_model(motif, starting_point_gs(max_iterations, tolerance, seqs, diagnostic_events, motif, f, geo_p, geo_delta, indicators));

//    cerr << "------------------GS------------------------" << endl;
//    f_liklihood = gibbs_sampling(max_iterations, tolerance, seqs, diagnostic_events, motif, 
//                     f, geo_p, geo_delta, indicators);

//      cout << print_model(motif, f, geo_p, geo_delta, regions, seqs, indicators); 
  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
