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
#include "MappedRead.hpp"
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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////
/////////////  A SET OF FUNCTIONS FOR READING THE MAPPING INFORMATION
/////////////                FOR DIAGNOSTIC EVENTS
/////////////

static bool
is_header_line(const string& line) {
  if (line.substr(0, 1) == "#")
    return true;
  static const char *browser_label = "browser";
  static const size_t browser_label_len = 7;
  for (size_t i = 0; i < browser_label_len; ++i)
    if (line[i] != browser_label[i])
      return false;
  return true;
}

static bool
is_track_line(const char *line) {
  static const char *track_label = "track";
  static const size_t track_label_len = 5;
  for (size_t i = 0; i < track_label_len; ++i)
    if (line[i] != track_label[i])
      return false;
  return true;
}

static vector<string>
split_extended_whitespace_quoted(string to_split) {

  static const char *non_word_chars = "\t\"'";
  to_split = smithlab::strip(to_split);
  std::vector<std::string> words;
  size_t start_pos = 0, end_pos = 0;
  const size_t length_of_to_split = to_split.length();
  while (start_pos < length_of_to_split) {
    /** find next position that is not a word character */
    end_pos = to_split.find_first_of(non_word_chars, end_pos);
    if (end_pos > to_split.length()) { /** If we hit the end: done */
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      break;
    }
    /** unescaped, unquoted white space: definitely a word delimiter */
    if (to_split[end_pos] == ' ' || to_split[end_pos] == '\t') {
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      end_pos = to_split.find_first_not_of(" \t", end_pos);
      start_pos = end_pos;
    }
    /** preserve whatever is being escaped; will become part of the
     *         current word */
    else if (to_split[end_pos] == '\\')
      end_pos = to_split.find_first_not_of(non_word_chars, end_pos + 2);
    else {
      const std::string current_delim = "\\" + to_split.substr(end_pos, 1);
      do { // slurp doubly- or singly-quoted string
        end_pos = to_split.find_first_of(current_delim, end_pos + 1);
        if (end_pos == std::string::npos) {
          end_pos = length_of_to_split;
          break;
        }
        if (to_split[end_pos] == '\\')
          ++end_pos;
        else
          break;
      } while (true);
      ++end_pos;
    }
    if (end_pos >= length_of_to_split) {
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      break;
    }
  }
  return words;
}

static size_t
convertString(const string& s) {
  istringstream buffer(s);
  size_t value;
  buffer >> value;
  return value;
}

static string
convertSizet(const size_t number) {
  stringstream ss;
  ss << number;
  return ss.str();
}

static void
convert_bowtie_extra(string &extra) {
  std::replace(extra.begin(), extra.end(), ',', ' ');
  extra.erase(std::remove(extra.begin(), extra.end(), ':'), extra.end());
}

static void
read_piranha_output(const string filename,
                    vector<GenomicRegion> &regions) {
  regions.clear();
  ReadBEDFile(filename, regions);
  for (size_t i = 0; i < regions.size(); ++i)
    regions[i].set_name("sequence_" + toa(i + 1));
}

static void
read_rmap_output(const vector<string> &parts,
                 vector<MappedRead> &regions) {
  if (parts.size() == 8) {
    string gr = parts[0] + "\t" + parts[1] + "\t" + parts[2] + "\t" + parts[3]
        + "\t" + parts[4] + "\t" + parts[5];
    string read = gr + "\t" + parts[6] + "\t2C>C";
    MappedRead mr(read.c_str());
    regions.push_back(mr);
  }
  else
    throw SMITHLABException("The mapped reads file is not correctly formatted!");
}

static void
read_bowtie_output(const vector<string> &parts,
                   vector<MappedRead> &regions) {
  if (parts.size() == 8) {
    vector<string> name(smithlab::split_whitespace_quoted(parts[0]));
    string extra = parts[7];
    convert_bowtie_extra(extra);
    string gr = parts[2] + "\t" + parts[3] + "\t"
        + convertSizet(convertString(parts[3]) + parts[4].length()) + "\t"
        + name[0] + "\t" + parts[6] + "\t" + parts[1];
    string read = gr + "\t" + parts[4] + "\t" + parts[7];
    MappedRead mr(read.c_str());
    regions.push_back(mr);
  }
}

static void
read_novoalign_output(const vector<string> &parts,
                      vector<MappedRead> &regions) {
  if (parts.size() == 14 && parts[4] == "U") {
    vector<string> name(smithlab::split_whitespace_quoted(parts[0]));
    string strand = ((parts[9] == "F") ? "+" : "-");
    string gr = parts[7].substr(1) + "\t" + parts[8] + "\t"
        + convertSizet(convertString(parts[8]) + parts[2].length()) + "\t"
        + name[0] + "\t" + parts[5] + "\t" + strand;
    string read = gr + "\t" + parts[2] + "\t" + parts[13];
    MappedRead mr(read.c_str());
    regions.push_back(mr);
  }
}

void
load_mapped_reads(const string &reads_file,
                  const string &mapper,
                  vector<MappedRead> &mapped_reads) {

  static const size_t buffer_size = 10000;

  std::ifstream in(reads_file.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + reads_file);

  std::ifstream::pos_type start_of_data = in.tellg();
  in.seekg(0, std::ios::end);
  std::ifstream::pos_type end_of_data = in.tellg();
  in.seekg(start_of_data);
  size_t current_done = 101;

  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw SMITHLABException("Line too long in file: " + reads_file);
    if (!is_header_line(buffer) && !is_track_line(buffer)) {
      vector<string> parts(split_extended_whitespace_quoted(buffer));
      if (mapper == "novoalign")
        read_novoalign_output(parts, mapped_reads);
      else if (mapper == "bowtie")
        read_bowtie_output(parts, mapped_reads);
      else if (mapper == "rmap")
        read_rmap_output(parts, mapped_reads);
      else
        throw SMITHLABException("The mapper was not recognized!");
    }
    size_t percent_done = static_cast<size_t>(in.tellg()) * 100 / end_of_data;
    if (percent_done != current_done) {
      std::cerr << "\r" << percent_done << "% completed..." << std::flush;
      current_done = percent_done;
    }
    in.peek();
  }
  in.close();
  std::cerr << std::endl << mapped_reads.size() << " uniquely mapped reads are loaded!" << std::endl;
}

////////////////////////////////////////////////////////////////////////

static void
expand_regions(vector<GenomicRegion> &regions,
               const size_t padding) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start(regions[i].get_start() - padding);
    regions[i].set_end(regions[i].get_end() + padding);
  }
}

static void
unexpand_regions(vector<GenomicRegion> &regions,
                 const size_t padding) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start(regions[i].get_start() + padding);
    regions[i].set_end(regions[i].get_end() - padding);
  }
}

static size_t
adjust_start_pos(const size_t orig_start,
                 const string &chrom_name) {
  static const double LINE_WIDTH = 50.0;
  // For the '>' and the '\n';
  const size_t name_offset = chrom_name.length() + 2;
  const size_t preceding_newlines = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  return orig_start + preceding_newlines + name_offset;
}

static size_t
adjust_region_size(const size_t orig_start,
                   const string &chrom_name,
                   const size_t orig_size) {
  static const double LINE_WIDTH = 50.0;
  const size_t preceding_newlines_start = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  const size_t preceding_newlines_end = static_cast<size_t>(std::floor(
      (orig_start + orig_size) / LINE_WIDTH));
  return (orig_size + (preceding_newlines_end - preceding_newlines_start));
}

static bool
test_chrom_file_format(const string filename) {
  static const double LINE_WIDTH = 50.0;
  static const double NEWLINES_TO_CHECK = 10.0;
  const size_t expected_lines =
  static_cast<size_t>(get_filesize(filename)/(LINE_WIDTH + 1.0)) + 1;
  

  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("could not open file: " + filename);

  string name;
  in >> name;
  if (name[0] != '>')
    throw SMITHLABException("bad format for chrom file: " + filename +\
                            " -- (expected name line with leading '>', " +\
                            "found " + name + " instead)");

  const size_t increment = expected_lines/NEWLINES_TO_CHECK;

  for (size_t i = 0; i < NEWLINES_TO_CHECK; ++i) {
    // get a position that should be a newline
    const size_t should_be_newline = 
      name.length() + i*increment*(LINE_WIDTH + 1);
    // go to that position
    in.seekg(should_be_newline, std::ios::beg);
    if (in.peek() != '\n' && in.peek() != '\r')
      throw SMITHLABException("bad format for chrom file: " + filename +\
                              " -- expected newline characters at 50nt " +\
                              "offsets; didn't find them.");
  }
  return true;
}

static void
extract_regions_chrom_fasta(const string &chrom_name,
                            const string &filename,
                            const vector<GenomicRegion> &regions,
                            vector<string> &sequences,
                            vector<string> &names) {
  test_chrom_file_format(filename);
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
    extract_regions_chrom_fasta(
        chrom_name, filenames[f_idx->second], regions[i], sequences, names);
  }
}

void
load_sequences(const string &chrom_dir,
               const size_t padding,
               const string &targets_file,
               vector<string> &names,
               vector<string> &sequences,
               vector<GenomicRegion> &targets) {

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
    if (!chrom_dir.empty()) {
      throw SMITHLABException("Provided chromosomes directory, but loading "
                              "sequences from fasta file; chromosomes "
                              "directory cannot be used.");
    }
    read_fasta_file(targets_file, names, sequences);
  } else {
    read_piranha_output(targets_file, targets);
    expand_regions(targets, padding);
    sort(targets.begin(), targets.end());
    if (chrom_dir.empty())
      throw SMITHLABException("Provided sequence genomic regions, but no "
                              "chromosomes directory to load the sequence "
                              "data from. Please specify a valid directory "
                              "containing the chromosome files.");
    else {
      if (!isdir(chrom_dir.c_str()))
        throw SMITHLABException(chrom_dir + " is not a valid directory");
    }
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
    throw SMITHLABException("cannot open input file " + structure_file);
  string s;
  while (getline(in, s)) {
    istringstream ss(s);
    vector<double> record;
    while (ss) {
      string s;
      if (!getline(ss, s, ','))
        break;
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
                    const string &outfile,
                    const size_t padding) {
  assert(!outfile.empty());
  std::ofstream out(outfile.c_str());
  for (size_t i = 0; i < sec_structure.size(); ++i) {
    assert(sec_structure.size() > 2 * padding);
    copy(
        sec_structure[i].begin() + padding, sec_structure[i].end() - padding,
        std::ostream_iterator<double>(out, ","));
    out << endl;
  }
}




