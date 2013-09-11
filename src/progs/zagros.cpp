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
using std::numeric_limits;



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////// 
///////////  CODE FOR FORMATTING MOTIF / MODEL OUTPUT BELOW HERE
///////////

static string 
format_site(const Model &model, const GenomicRegion &region, 
	    const string &seq, const size_t site_pos) {
  const size_t width = model.get_model_size();
  std::ostringstream ss;
  ss << "BS\t" << seq.substr(site_pos, width) << "; "
     << assemble_region_name(region) << "; " 
     << site_pos << "; "
     << width << ";  ;"
     << (region.pos_strand() ? "p;" : "n;");
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
format_motif(const Model &model, const string &motif_name,
	     const vector<GenomicRegion> &targets, 
	     const vector<string> &sequences,
	     const vector<vector<double> > &indicators, 
	     const vector<double> &zoops_i) {
  
  assert(sequences.size() == indicators.size() &&
	 zoops_i.size() == sequences.size());
  
  std::ostringstream ss;
  ss << format_motif_header(motif_name) << endl;
  
  vector<vector<double> > tmp_m = model.getM();
  for (size_t n = 0; n < sequences.size(); n++) {
    double max_X = -1;
    int max_i = -1;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < model.get_model_size(); ++j)
      tmp_m[j][model.base2int_RNA(sequences[n][max_i + j])] += 1;
  }
  
  for (size_t j = 0; j < tmp_m.size(); j++) {
    // AS: this is crazy below and will not work for motifs of width 11
    ss << "0" << j + 1;
    for (size_t b = 0; b < smithlab::alphabet_size; ++b)
      ss << '\t' << static_cast<int>(tmp_m[j][b]);
    ss << endl;
  }
  ss << "XX" << endl
     << "AT\tGEO_P=" << model.getP() << endl
     << "XX" << endl;
  
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
      if (zoops_i[n] > model.zoops_threshold)
	ss << format_site(model, targets[n], sequences[n], site_pos) << "\t" << zoops_i[n] << endl;
  }
  ss << "XX" << endl
     << "//" << endl;
  
  return ss.str();
}



int main(int argc, const char **argv) {

  try {

    static const size_t flanking_regions_size = 20;
    
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
    OptionParser opt_parse(strip_path(argv[0]), "", 
			   "<target regions/sequences>");
    opt_parse.add_opt("output", 'o', "name of output directory "
		      "(default: zagros_output)", false, outdir);
    opt_parse.add_opt("width", 'w', "motif width", false, motif_width);
    opt_parse.add_opt("chrom", 'c', "directory with chrom files (FASTA format)",
		      false, chrom_dir);
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
    
    mkdir(outdir.c_str(), 0750);
    
    string dirname, base_name, suffix;
    parse_dir_baseanme_suffix(targets_file, dirname, base_name, suffix);
    const string base_file(path_join(outdir, base_name));
    
    const size_t padding =
      structure_information_file.empty() ? 0 : flanking_regions_size;
    
    if (!structure_information_file.empty()) {
      std::ifstream in(structure_information_file.c_str(), std::ios::binary);
      if (!in)
        throw SMITHLABException("cannot open input file " + 
				structure_information_file);
    }
    
    //Vectors to store primary information from the data
    vector<string> seqs, names;
    vector<GenomicRegion> targets;
    vector<vector<size_t> > diagnostic_events;
    vector<double> secondary_structure;
    load_sequences(chrom_dir, padding, targets_file, names, seqs, targets);
    
    vector<double> has_motif(seqs.size(), 0.5);
    vector<vector<double> > indicators;
    for (size_t i = 0; i < seqs.size(); ++i) {
      const size_t window_size =
	targets.empty() ? seqs[i].length() : targets[i].get_width();
      const size_t n_pos = window_size - motif_width + 1;
      indicators.push_back(vector<double>(n_pos, 1.0/n_pos));
    }
    
    Model model(motif_width);
    model.expectation_maximization(seqs, diagnostic_events, 
				   secondary_structure, indicators, 
				   has_motif, structure_information_file, 
				   use_de_information, base_file,
				   max_iterations, tolerance);
    
    const string output_model(base_file + ".mat");
    std::ofstream out(output_model.c_str());
    out << format_motif(model, "ME_EM", targets, seqs, 
			indicators, has_motif) << endl;
    
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
