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
#include "ExtendedGenomicRegion.hpp"

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
    bool use_sequence_information = false;
    bool use_structure_information = false;
    bool use_de_information = false;
    size_t max_de = std::numeric_limits<size_t>::max();
    size_t min_cluster_size = 1;

    string reads_file;

    const size_t max_iterations = 10;
    const double tolerance = 1e-10;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<targets>");
    opt_parse.add_opt(
        "output", 'o', "Name of output directory (default: zagros_output)",
        false, outdir);
    opt_parse.add_opt("width", 'w', "motif width", false, motif_width);
    opt_parse.add_opt(
        "chrom", 'c', "directory with chrom files (FASTA format)", true,
        chrom_dir);
    opt_parse.add_opt(
        "experiment",
        'e',
        "The type of experiment: iCLIP, hCLIP, pCLIP (default: no diagnostic events considered)",
        false, experiment);
    opt_parse.add_opt("mapper", 'm', "Mapped reads format: novoalign, "
        "bowtie, rmap, piranha (default: none)", false, mapper);
    opt_parse.add_opt(
        "sequence", 's', "Use the sequence information", false,
        use_sequence_information);
    opt_parse.add_opt("reads", 'r', "Mapped reads file", false, reads_file);
    opt_parse.add_opt(
        "structure", 't', "Use the structure information", false,
        use_structure_information);
    opt_parse.add_opt(
        "diagnostic_events", 'd', "Use the diagnostic events information",
        false, use_de_information);
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
    parse_dir_baseanme_suffix(targets_file,dirname,
        base_name, suffix);
    //Create the base file name for output files
    string base_file = outdir + base_name;

    //Vectors to store primary information from the data
    vector<string> seqs;
    vector<string> names;
    vector<GenomicRegion> de_regions;
    vector<GenomicRegion> targets;
    vector<vector<size_t> > diagnostic_events;

    //Reading the targets
    cerr << "Reading target regions...";
    IO::read_piranha_output(targets_file, targets);
    cerr << "done!" << endl;
    if (targets.size() == 0)
      throw SMITHLABException("No targets found...");
    //Sorting the target file and storing the file in output directory
    sort(targets.begin(), targets.end(), region_less());

    //Extracting the target sequences
    if (use_structure_information) {
      IO::expand_regions(targets);
      IO::extract_regions_fasta(chrom_dir, targets, seqs, names);
      IO::unexpand_regions(targets);
    } else
      IO::extract_regions_fasta(chrom_dir, targets, seqs, names);

    unordered_map<string, size_t> names_table;
    IO::make_sequence_names(names, seqs, targets, names_table);

    //Reading the mapped reads
    if (!reads_file.empty()) {
      vector<ExtendedGenomicRegion> mapped_reads;

      cerr << "Reading mapped reads..." << endl;
      if (mapper == "rmap")
        IO::read_rmap_output(reads_file, mapped_reads);
      else if (mapper == "novoalign")
        IO::read_novoalign_output(reads_file, mapped_reads);
      else if (mapper == "bowtie")
        IO::read_bowtie_output(reads_file, mapped_reads);
      else
        throw SMITHLABException("Mapper not recognized, please choose from "
            "the options (rmap, novoalign or bowtie)");

      //Making the regions and extracting the diagnostic events from the mapped
      //reads
      cerr << "Processing mapped reads..." << endl;
      vector<GenomicRegion> regions;
      IO::make_inputs(
          mapped_reads, regions, de_regions, experiment, max_de,
          min_cluster_size);

      if (regions.size() == 0)
        throw SMITHLABException("No reads found...");
      //Sorting the diagnostic events file and storing the file in output directory
      sort(de_regions.begin(), de_regions.end(), region_less());

      IO::sift(targets, de_regions);
      IO::load_diagnostic_events(
          de_regions, names_table, targets, diagnostic_events);
    }

    vector<double> has_motif(targets.size(), 1.0);
    vector<vector<double> > indicators;
    for (size_t i = 0; i < targets.size(); ++i) {
      const size_t n_pos = targets[i].get_width() - motif_width + 1;
      indicators.push_back(vector<double>(n_pos, 1.0 / n_pos));
    }

    Model model(motif_width);

    cerr << "Fitting started..." << endl;
    model.expectation_maximization(
        max_iterations, tolerance, targets, seqs, diagnostic_events, indicators,
        has_motif, use_sequence_information, use_structure_information,
        use_de_information, base_file);

    cerr << "Preparing output...";
    IO::save_input_files(seqs, targets, de_regions, base_file);
    if (use_de_information)
      model.prepare_output(seqs, indicators, diagnostic_events, base_file);
    else
      model.prepare_output(seqs, indicators, base_file);
    string output_model = base_file + ".mat";
    ofstream outf(output_model.c_str());
    outf << model.print_model("DE_EM", targets, seqs, indicators);
    cerr << "done!" << endl;
  } catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
