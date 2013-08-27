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

using smithlab::alphabet_size;

int main(int argc, const char **argv) {

	try {

		bool VERBOSE = false;
		string outdir = "zagros_output";
		string chrom_dir;
		size_t motif_width = 8;
		string experiment = "none";
		string mapper = "piranha";
		bool use_sequence_information = false;
		bool use_structure_information = false;
		size_t max_de = 0;
		size_t min_cluster_size = 100;

		const size_t max_iterations = 10;
		const double tolerance = 1e-10;

		/****************** COMMAND LINE OPTIONS ********************/
		OptionParser opt_parse(strip_path(argv[0]), "", "<mapped reads/peaks>");
		opt_parse.add_opt("output", 'o',
				"Name of output directory (default: zagros_output)", false,
				outdir);
		opt_parse.add_opt("width", 'w', "motif width", false, motif_width);
		opt_parse.add_opt("chrom", 'c',
				"directory with chrom files (FASTA format)", true, chrom_dir);
		opt_parse.add_opt("experiment", 'e',
				"The type of experiment: iCLIP, hCLIP, pCLIP (default: no diagnostic events considered)",
				false, experiment);
		opt_parse.add_opt("mapper", 'm',
				"Mapped reads format: novoalign, bowtie, rmap, piranha (default: piranha)",
				false, mapper);
		opt_parse.add_opt("sequence", 's', "Use the sequence information",
				false, use_sequence_information);
		opt_parse.add_opt("structure", 't', "Use the structure information",
				false, use_structure_information);
		opt_parse.add_opt("diagnostic_events", 'd',
				"Maximum number of diagnostic events per sequence", false,
				max_de);
		opt_parse.add_opt("verbose", 'v', "print more run info", false,
				VERBOSE);
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
			cerr << opt_parse.help_message() << endl;
			return EXIT_SUCCESS;
		}

		const string reads_file(leftover_args.front());

		/****************** END COMMAND LINE OPTIONS *****************/

		vector<string> seqs;
		vector<string> names;
		vector<GenomicRegion> regions;
		vector<GenomicRegion> de_regions;

		vector<ExtendedGenomicRegion> mapped_reads;
		if (mapper == "rmap" || mapper == "novoalign" || mapper == "bowtie")
			cerr << "Reading mapped reads..." << endl;
		else
			cerr << "Reading input regions..." << endl;

		if (mapper == "rmap")
			read_rmap_output(reads_file, mapped_reads);
		else if (mapper == "novoalign")
			read_novoalign_output(reads_file, mapped_reads);
		else if (mapper == "bowtie")
			read_bowtie_output(reads_file, mapped_reads);
		else {
			read_piranha_output(reads_file, mapped_reads);
			min_cluster_size = 0;
		}

		cerr << "Processing mapped reads..." << endl;
		IO::make_inputs(mapped_reads, regions, de_regions, experiment, max_de,
				min_cluster_size);

		if (outdir.find_last_of("/\\") != outdir.length() - 1)
			outdir += "/";
		mkdir(outdir.c_str(), 0750);

		size_t base_index = reads_file.find_last_of("/\\");
		size_t suffix_index = reads_file.find_last_of(".");
		if (suffix_index <= base_index)
			suffix_index = reads_file.length() - 1;
		string base_file = outdir
				+ reads_file.substr(base_index + 1,
						suffix_index - base_index - 1);

		if (regions.size() == 0)
			throw BEDFileException("No significant clusters found...");

		string regions_outfile = base_file + ".bed";
		IO::sort_regions(regions, regions_outfile);
		regions.clear();
		ReadBEDFile(regions_outfile, regions);

                IO::expand_regions(regions);
		IO::extract_regions_fasta(chrom_dir, regions, seqs, names);
		IO::unexpand_regions(regions);

		unordered_map<string, size_t> names_table;
		IO::make_sequence_names(names, seqs, regions, names_table);

		vector<vector<size_t> > diagnostic_events;
		IO::load_diagnostic_events(de_regions, names_table, regions,
				diagnostic_events);

		vector<vector<double> > indicators;
		for (size_t i = 0; i < regions.size(); ++i) {
			const size_t n_pos = regions[i].get_width() - motif_width + 1;
			indicators.push_back(vector<double>(n_pos, 1.0 / n_pos));
		}

		Model model(motif_width, regions);
		model.set_delta(0);
//		model.find_delta(seqs, regions, diagnostic_events);
//		cout << model.get_delta() << endl;

    IO::save_input_files(seqs, regions, de_regions, base_file);

    model.expectation_maximization(max_iterations, tolerance, regions, seqs,
				diagnostic_events, indicators, use_sequence_information,
				use_structure_information, max_de, base_file);

		model.prepare_output(seqs, indicators, diagnostic_events, base_file);
		string output_model = base_file + ".mat";
		ofstream outf(output_model.c_str());
		outf << model.print_model("DE_EM", regions, seqs, indicators);
	} catch (const SMITHLABException &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	} catch (std::bad_alloc &ba) {
		cerr << "ERROR: could not allocate memory" << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
