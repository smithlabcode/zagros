/*
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith, Emad Bahrami-Samani, Philip J. Uren
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

// stl includes
#include <string>
#include <vector>

// smithlab common code includes
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

// bring these into the default namespace
using std::cerr;
using std::endl;
using std::string;
using std::vector;

/**
 * \brief The main entry point for the Scan program.
 * \param argc count of the number of arguments passed to the program.
 * \param argv array of char*, each one is an argument passed to the program.
 */
int main(int argc, const char **argv) {
  const string about = "This program scans sequences for occurrences of a "
                       "motif discovered by Zagros. For each input "
                       "sequence, it outputs the sequence and (relative) "
                       "position of the best match to the input motif.";

  string infn = "";     // this is the motif file output by Zagros.
  string seqfn = "";    // these are the sequence to scan, in fasta format.
  string strfn = "";    // these are the base-pair-prob. scores for the seqs
  string desfn = "";    // these are the diagnostic event counts for the seqs
  string outfn = "";

  bool VERBOSE = false;

  try {
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), about,
                           "<target_regions/sequences>");
    opt_parse.add_opt("sequences", 's', "sequences to scan (fasta format)",
                      OptionParser::REQUIRED, seqfn);
    opt_parse.add_opt("structure", 't', "structure information file",
                      OptionParser::OPTIONAL, strfn);
    opt_parse.add_opt("diagnostic_events", 'd', "diagnostic events "
                                                "information file",
                      OptionParser::OPTIONAL, desfn);
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      OptionParser::OPTIONAL, outfn);
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
    if (leftover_args.size() != 1) {
      cerr << "This program requires exactly one input file (the motif "
           << "file, as output by Zagros), found "
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

    cerr << "SORRY, THIS PROGRAM IS NOT IMPLEMENTED YET." << endl;
  }
  catch (const SMITHLABException &e) {
    cerr << "ERROR: " << e.what();
    cerr << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
}

