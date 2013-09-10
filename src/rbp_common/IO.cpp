/*    sortbed: a program for sorting BED format files
 *    Copyright (C) 2008 University of Southern California and
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
#include "GenomicRegion.hpp"
#include "Model.hpp"
#include "IO.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::accumulate;
using std::ostream;
using std::pair;
using std::make_pair;
using std::sort;
using std::numeric_limits;

using smithlab::alphabet_size;

using std::tr1::unordered_map;

size_t convertString(const string& s) {
  istringstream buffer(s);
  size_t value;
  buffer >> value;
  return value;
}

string convertSizet(const size_t number) {
  stringstream ss;
  ss << number;
  return ss.str();
}

IO::IO() {
}

void IO::load_sequences(vector<string> &names, vector<string> &sequences,
    vector<GenomicRegion> &targets, const string chrom_dir,
    const size_t padding, const string targets_file) {

  std::ifstream in(targets_file.c_str(), std::ios::binary);
  if (!in) {
    throw SMITHLABException("cannot open input file " + string(targets_file));
  }

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  char buffer[INPUT_BUFFER_SIZE + 1];
  in.getline(buffer, INPUT_BUFFER_SIZE);
  char first_char = buffer[0];
  in.close();
  if (first_char == '>') {
    if (padding != 0)
      throw SMITHLABException(
          "Input the genomic regions, if you wish to use the secondary structure!");
    read_fasta_file(targets_file, names, sequences);
  } else {
    read_piranha_output(targets_file, targets);
    expand_regions(targets, padding);
    sort(targets.begin(), targets.end(), region_less());
    if (chrom_dir == "")
      throw SMITHLABException(
          "Input a valid directory containing the chromosome files!");
    extract_regions_fasta(chrom_dir, targets, sequences, names);
    unexpand_regions(targets, padding);
  }
}

void IO::read_piranha_output(const string filename,
    vector<GenomicRegion> &regions) {

  vector<GenomicRegion> the_regions;
  ReadBEDFile(filename, the_regions);
  for (size_t i = 0; i < the_regions.size(); i += 1) {
    string name = "sequence_" + convertSizet(i + 1);
    GenomicRegion gr(
        the_regions[i].get_chrom(), the_regions[i].get_start(),
        the_regions[i].get_end(), name, the_regions[i].get_score(),
        the_regions[i].get_strand());
    regions.push_back(gr);
  }
}

void IO::extract_regions_fasta(const string &dirname,
    const vector<GenomicRegion> &regions_in, vector<string> &sequences,
    vector<string> &names) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<GenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
        chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw SMITHLABException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(
        chrom_name, filenames[f_idx->second], regions[i], sequences, names);
  }
}

void IO::extract_regions_chrom_fasta(const string &chrom_name,
    const string &filename, const vector<GenomicRegion> &regions,
    vector<string> &sequences, vector<string> &names) {

  std::ifstream in(filename.c_str());
  for (vector<GenomicRegion>::const_iterator i(regions.begin());
      i != regions.end(); ++i) {

    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size = adjust_region_size(
        orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(
        buffer, buffer + region_size,
        std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    names.push_back(i->get_name());
    std::transform(
        sequences.back().begin(), sequences.back().end(),
        sequences.back().begin(), std::ptr_fun(&toupper));
    if (i->neg_strand())
      revcomp_inplace(sequences.back());
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}

size_t IO::adjust_start_pos(const size_t orig_start, const string &chrom_name) {
  static const double LINE_WIDTH = 50.0;
  const size_t name_offset = chrom_name.length() + 2; // For the '>' and the '\n';
  const size_t preceding_newlines = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  return orig_start + preceding_newlines + name_offset;
}

size_t IO::adjust_region_size(const size_t orig_start, const string &chrom_name,
    const size_t orig_size) {
  static const double LINE_WIDTH = 50.0;
  const size_t preceding_newlines_start = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  const size_t preceding_newlines_end = static_cast<size_t>(std::floor(
      (orig_start + orig_size) / LINE_WIDTH));
  return (orig_size + (preceding_newlines_end - preceding_newlines_start));
}

void IO::expand_regions(vector<GenomicRegion> &regions, const size_t padding) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start(regions[i].get_start() - padding);
    regions[i].set_end(regions[i].get_end() + padding);
  }
}

void IO::unexpand_regions(vector<GenomicRegion> &regions,
    const size_t padding) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start(regions[i].get_start() + padding);
    regions[i].set_end(regions[i].get_end() - padding);
  }
}

std::string IO::print_model(const Model &model, const string &motif_name,
    const vector<GenomicRegion> &targets, const vector<string> &sequences,
    const vector<vector<double> > &indicators) {
  stringstream ss;

  const size_t N = sequences.size();
  if (N <= 0) {
    stringstream ss;
    ss << "Building motif alignment failed. Reason: sequences vector is empty";
    throw SMITHLABException(ss.str());
  }

  if (N != indicators.size()) {
    stringstream ss;
    ss << "Building motif alignment failed. Reason: expected " << N
        << "indicator vectors, got only " << indicators.size();
    throw SMITHLABException(ss.str());
  }

  ss << "AC\t" << motif_name << endl;
  ss << "XX" << endl;
  ss << "TY\tMotif" << endl;
  ss << "XX" << endl;
  ss << "P0\tA\tC\tG\tT" << endl;

  double max_X = -1;
  int max_i = -1;

  vector<vector<double> > tmp_m = model.getM();
  for (size_t n = 0; n < sequences.size(); n++) {
    max_X = -1;
    max_i = -1;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < model.get_model_size(); ++j) {
      tmp_m[j][model.base2int_RNA(sequences[n][max_i + j])] += 1;
    }
  }

  for (size_t j = 0; j < tmp_m.size(); j++) {
    ss << "0" << j + 1 << "\t";
    for (size_t b = 0; b < alphabet_size - 1; b++)
      ss << (int) (tmp_m[j][b]) << "\t";
    ss << (int) (tmp_m[j][alphabet_size - 1]) << endl;
  }
  ss << "XX" << endl;
  ss << "AT\tGEO_P=" << model.getP() << endl;
  ss << "XX" << endl;

  for (size_t n = 0; n < N; n++) {
    max_X = -1;
    max_i = 0;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    if (!targets.empty())
      ss << "BS\t" << sequences[n].substr(max_i, model.get_model_size()) << "; "
          << targets[n].get_chrom() << ":" << targets[n].get_start() << "-"
          << targets[n].get_end() << "; " << max_i + 1 << "; "
          << model.get_model_size() << ";  ;"
          << ((targets[n].get_strand() == '+') ? "p; " : "n;") << endl;
  }
  ss << "XX" << endl;
  ss << "//" << endl;

  return ss.str();
}

