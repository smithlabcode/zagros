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

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <numeric>
#include <tr1/unordered_map>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "Model.hpp"
#include "IO.hpp"
#include "RNG.hpp"

using std::string;
using std::vector;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ostream;
using std::sort;

using smithlab::alphabet_size;

using std::tr1::unordered_map;

static void
read_piranha_output(const string filename,
                    vector<GenomicRegion> &regions) {
  regions.clear();
  ReadBEDFile(filename, regions);
  for (size_t i = 0; i < regions.size(); ++i)
    regions[i].set_name("sequence_" + toa(i + 1));
}

static void
expand_regions(vector<GenomicRegion> &regions, const size_t padding) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start(regions[i].get_start() - padding);
    regions[i].set_end(regions[i].get_end() + padding);
  }
}

static void
unexpand_regions(vector<GenomicRegion> &regions, const size_t padding) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start(regions[i].get_start() + padding);
    regions[i].set_end(regions[i].get_end() - padding);
  }
}

static size_t
adjust_start_pos(const size_t orig_start, const string &chrom_name) {
  static const double LINE_WIDTH = 50.0;
  // For the '>' and the '\n';
  const size_t name_offset = chrom_name.length() + 2;
  const size_t preceding_newlines =
    static_cast<size_t>(std::floor(orig_start/LINE_WIDTH));
  return orig_start + preceding_newlines + name_offset;
}

static size_t
adjust_region_size(const size_t orig_start, const string &chrom_name,
                   const size_t orig_size) {
  static const double LINE_WIDTH = 50.0;
  const size_t preceding_newlines_start =
    static_cast<size_t>(std::floor(orig_start/LINE_WIDTH));
  const size_t preceding_newlines_end =
    static_cast<size_t>(std::floor((orig_start + orig_size)/LINE_WIDTH));
  return (orig_size + (preceding_newlines_end - preceding_newlines_start));
}

static void
extract_regions_chrom_fasta(const string &chrom_name,
                            const string &filename,
                            const vector<GenomicRegion> &regions,
                            vector<string> &sequences,
                            vector<string> &names) {

  std::ifstream in(filename.c_str());
  for (vector<GenomicRegion>::const_iterator i(regions.begin());
       i != regions.end(); ++i) {

    const size_t orig_start_pos = i->get_start();
    const size_t orig_region_size = i->get_end() - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size = adjust_region_size(
        orig_start_pos, chrom_name, orig_region_size);
    // It's a "size_t", so of course it will be >= 0!!!
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(buffer, buffer + region_size,
                   std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    names.push_back(i->get_name());
    std::transform(sequences.back().begin(), sequences.back().end(),
                   sequences.back().begin(), std::ptr_fun(&toupper));
    if (i->neg_strand())
      revcomp_inplace(sequences.back());
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}

static void
extract_regions_fasta(const string &dirname,
                      const vector<GenomicRegion> &regions_in,
                      vector<string> &sequences,
                      vector<string> &names) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<GenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    unordered_map<string, size_t>::const_iterator f_idx =
        chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw SMITHLABException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(chrom_name, filenames[f_idx->second],
                                regions[i], sequences, names);
  }
}

void
load_sequences(const string &chrom_dir, const size_t padding,
               const string &targets_file, vector<string> &names,
               vector<string> &sequences, vector<GenomicRegion> &targets) {

  std::ifstream in(targets_file.c_str(), std::ios::binary);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(targets_file));

  string buffer;
  getline(in, buffer);
  in.close();

  if (buffer[0] == '>') {
    if (padding != 0)
      throw SMITHLABException("Input the genomic regions, "
                              "if you wish to use the "
                              "secondary structure!");
    read_fasta_file(targets_file, names, sequences);
  } else {
    read_piranha_output(targets_file, targets);
    expand_regions(targets, padding);
    sort(targets.begin(), targets.end());
    if (chrom_dir.empty())
      throw SMITHLABException("Input a valid directory containing "
          "the chromosome files!");
    extract_regions_fasta(chrom_dir, targets, sequences, names);
    unexpand_regions(targets, padding);
  }
}

bool
seq_and_structure_are_consistent(const vector<string> &seqs,
                                 const vector<vector<double> > &sec_structure) {
  if (seqs.size() != sec_structure.size())
    return false;
  for (size_t i = 0; i < seqs.size(); ++i)
    if (seqs[i].length() != sec_structure[i].size())
      return false;
  return true;
}

void
load_structures(const string structure_file,
                vector<vector<double> > &structures) {

  std::ifstream in(structure_file.c_str(), std::ios::binary);
  if (!in)
    throw SMITHLABException("cannot open input file " +
                            structure_file);
  string s;
  while (getline(in, s)) {
    istringstream ss(s);
    vector<double> record;
    while (ss) {
      string s;
      if (!getline(ss, s, ',')) break;
      double d;
      stringstream s2d(s);
      s2d >> d;
      record.push_back(d);
    }
    structures.push_back(record);
  }
}

void
save_structure_file(const vector<vector<double> > &sec_structure,
                    const string &outfile, const size_t padding) {
  assert(!outfile.empty());
  std::ofstream out(outfile.c_str());
  for (size_t i = 0; i < sec_structure.size(); ++i) {
    assert(sec_structure.size() > 2*padding);
    copy(sec_structure[i].begin() + padding, sec_structure[i].end() - padding,
         std::ostream_iterator<double>(out, ","));
    out << endl;
  }
}
