/*    Functions for Zagros I/O
 *
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

/******************************************************************************
 *                STATIC HELPER FUNCTIONS FOR IO OPERATIONS
 ******************************************************************************/

/**
 * \brief convert a string to a size_t
 * \param s the string to convert
 */
static size_t
convertString(const string& s) {
  istringstream buffer(s);
  size_t value;
  buffer >> value;
  return value;
}

/**
 * \brief convert a size_t to a string
 * \param number the size_t to convert
 */
static string
convertSizet(const size_t number) {
  stringstream ss;
  ss << number;
  return ss.str();
}

/**
 *  * \brief convert a string to a size_t
 *   * \param   s   the string to convert
 *    * \return      the size_t corresponding to s
 *     */
static size_t
stringToSize_t(const string& s) {
  std::istringstream buffer(s);
  size_t value;
  buffer >> value;
  return value;
}

/**
 * \brief TODO
 * \param line  TODO
 * \return TODO
 */
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

/**
 * \brief TODO
 * \param line  TODO
 * \return TODO
 */
static bool
is_track_line(const char *line) {
  static const char *track_label = "track";
  static const size_t track_label_len = 5;
  for (size_t i = 0; i < track_label_len; ++i)
    if (line[i] != track_label[i])
      return false;
  return true;
}

/**
 * \brief TODO
 * \param to_split TODO
 * \return TODO
 * TODO -- PJU: This looks identical to the function in smithlab_utils.
 *              Why is it duplicated here?
 */
static vector<string>
split_extended_whitespace_quoted(string to_split) {

  static const char *non_word_chars = "\t";
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



static void
convert_bowtie_extra(string &extra, size_t read_length, char strand) {
  std::replace(extra.begin(), extra.end(), ',', ' ');
  extra.erase(std::remove(extra.begin(), extra.end(), ':'), extra.end());
  if (strand == '-') {
    vector<string> parts(smithlab::split_whitespace_quoted(extra));
    string result = "";
    size_t position;
    string refBase;
    string readBase;
    for (size_t i = 0; i < parts.size(); ++i) {
      size_t idx = parts[i].find(">");
      if (idx != string::npos) {
        position = stringToSize_t(parts[i].substr(0, idx - 1));
        refBase = parts[i].substr(idx - 1, 1);
        readBase = parts[i].substr(idx + 1, 1);
        result = result + " " + convertSizet(read_length - position) + refBase + ">" + readBase;
      }
      else {
        idx = parts[i].find("+");
        if (idx != string::npos) {
          position = stringToSize_t(parts[i].substr(0, idx));
          readBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
          refBase = "-";
          result = result + " " + convertSizet(read_length - position) + "+" + readBase;
        } else {
          idx = parts[i].find("-");
          if (idx != string::npos) {
            position = stringToSize_t(parts[i].substr(0, idx));
            refBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
            readBase = "-";
            result = result + " " + convertSizet(read_length - position) + "-" + refBase;
          }
        }
      }
    }
    extra = result.substr(1);
  } else
  if (strand == '+') {
    vector<string> parts(smithlab::split_whitespace_quoted(extra));
    string result = "";
    size_t position;
    string refBase;
    string readBase;
    for (size_t i = 0; i < parts.size(); ++i) {
      size_t idx = parts[i].find(">");
      if (idx != string::npos) {
        position = stringToSize_t(parts[i].substr(0, idx - 1));
        refBase = parts[i].substr(idx - 1, 1);
        readBase = parts[i].substr(idx + 1, 1);
        result = result + " " + convertSizet(position + 1) + refBase + ">" + readBase;
      }
      else {
        idx = parts[i].find("+");
        if (idx != string::npos) {
          position = stringToSize_t(parts[i].substr(0, idx));
          readBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
          refBase = "-";
          result = result + " " + convertSizet(position + 1) + "+" + readBase;
        } else {
          idx = parts[i].find("-");
          if (idx != string::npos) {
            position = stringToSize_t(parts[i].substr(0, idx));
            refBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
            readBase = "-";
            result = result + " " + convertSizet(position + 1) + "-" + refBase;
          }
        }
      }
    }
    extra = result.substr(1);
  }
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
    string read = gr + "\t" + parts[6] + "\t" + parts[7];
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
    char strand = ((parts[1] == "+") ? '+' : '-');
    string extra = parts[7];
    convert_bowtie_extra(extra, parts[4].length(), strand);
    string gr = parts[2] + "\t" + parts[3] + "\t"
        + convertSizet(convertString(parts[3]) + parts[4].length()) + "\t"
        + name[0] + "\t" + parts[6] + "\t" + parts[1];
    string read = gr + "\t" + parts[4] + "\t" + parts[7];
    MappedRead mr(read.c_str());
    regions.push_back(mr);
  }
}

/******************************************************************************
 *             STATIC HELPER FUNCTIONS FOR PARSING MAPPED READS
 ******************************************************************************/
/**
 * \brief determine whether a tokenized novoalign read is uniquely mapping
 */
static bool
isUniqueMapper_novoAlign(const vector<string> &parts) {
  // TODO fix magic numbers
  return (parts.size() >= 13 && parts[4] == "U");
}

/**
 * \brief determine whether a tokenized read from a given short read mapper is
 *        uniquely mapping.
 */
static bool
isUniqueMapper(const vector<string> &tokens, const string &mapper) {
  if (mapper == "novoalign") return isUniqueMapper_novoAlign(tokens);
  else throw SMITHLABException("The mapper '" + mapper +\
                               "' was not recognized!");
}

/**
 * \brief parse a novoalign read string that has already been tokenized and
 *        populate a MappedRead object. The read string must represent a
 *        uniquely mapping read. Note that, rather than putting the quality
 *        scores into the mr.scr field, we are putting the novoalign mismatch
 *        string into this field.
 * \param read TODO
 * \param mr   TODO
 * \throw SMITHLABException if the read string is not in the expected format.
 */
static void
parseMappedRead_novoAlign(const vector<string> &parts, MappedRead &mr) {
  // TODO fix magic numbers 13 and 4
  if (parts.size() >= 13 && parts[4] == "U") {
    vector<string> name(smithlab::split_whitespace_quoted(parts[0]));
    char strand = ((parts[9] == "F") ? '+' : '-');
    mr.r.set_chrom(parts[7].substr(1));
    mr.r.set_start(convertString(parts[8]));
    mr.r.set_end(convertString(parts[8]) + parts[2].length());
    mr.r.set_name(name[0]);
    mr.r.set_score(convertString(parts[5]));
    mr.r.set_strand(strand);
    mr.seq = parts[2];
    if (parts.size() >= 14) mr.scr = parts[13];
  } else {
    throw SMITHLABException("failed to parse novoalign read. Unexpected "
                            "format. Can only parse uniquely mapping reads");
  }
}

static void
parseMappedRead_bowtie(const vector<string> &parts, MappedRead &mr) {
  if (parts.size() >= 7) {
    vector<string> name(smithlab::split_whitespace_quoted(parts[0]));
    char strand = ((parts[1] == "+") ? '+' : '-');
    mr.r.set_chrom(parts[2]);
    mr.r.set_start(convertString(parts[3]) + 1);
    mr.r.set_end(convertString(parts[3]) + 1 + parts[4].length());
    mr.r.set_name(name[0]);
    mr.r.set_score(convertString(parts[6]));
    mr.r.set_strand(strand);
    mr.seq = parts[4];
    string extra;
    if (parts.size() == 8) {
      extra = parts[7];
      convert_bowtie_extra(extra, parts[4].length(), strand);
    }
    if (parts.size() == 8) mr.scr = extra;
  } else {
    throw SMITHLABException("failed to parse novoalign read. Unexpected "
                            "format. Can only parse uniquely mapping reads");
  }
}

static void
parseMappedRead_rmap(const vector<string> &parts, MappedRead &mr) {
  if (parts.size() >= 6) {
    char strand = ((parts[5] == "+") ? '+' : '-');
    mr.r.set_chrom(parts[0]);
    mr.r.set_start(convertString(parts[1]));
    mr.r.set_end(convertString(parts[2]));
    mr.r.set_name(parts[3]);
    mr.r.set_score(convertString(parts[4]));
    mr.r.set_strand(strand);
    if (parts.size() == 8) {
      mr.seq = parts[6];
      mr.scr = parts[7];
    } else {
      mr.seq = "NA";
      mr.scr = "NA";
    }
  } else {
    throw SMITHLABException("failed to parse novoalign read. Unexpected "
                            "format. Can only parse uniquely mapping reads");
  }
}


/**
 * \brief parse a tokenized novoalign read-string, create a mapped read from
 *        it and add the new MR to the provided vector. The tokenized read must
 *        be uniquely mapping.
 * \param ps the tokenized read ps[i] is the ith token from the (tab delimited)
 *           string representation.
 * \param rs the vector of mapped reads to add the resultant MappedRead to.
 */
static void
parseMappedRead_novoAlign(const vector<string> &ps, vector<MappedRead> &rs) {
  MappedRead mr;
  parseMappedRead_novoAlign(ps, mr);
  rs.push_back(mr);
}

/**
 * \brief parse a tokenized mapped read and populate the provide MappedRead
 *        object.
 * \param parts         TODO
 * \param mapper        TODO
 * \param mr            TODO
 */
static void
parseMappedRead(const vector<string> &parts, const string &mapper,
                MappedRead &mr) {
  if (mapper == "novoalign") parseMappedRead_novoAlign(parts, mr);
  else if (mapper == "bowtie") parseMappedRead_bowtie(parts, mr);
  else if (mapper == "rmap") parseMappedRead_rmap(parts, mr);
  else throw SMITHLABException("The mapper '" + mapper +\
                                                     "' was not recognized!");
}

/**
 * \brief parse a tokenized mapped read and add the resulting MappedRead object
 *        to a vector of mapped reads.
 * \param parts          TODO
 * \param mapper        TODO
 * \param mapped_reads  TODO
 */
static void
parseMappedRead(const vector<string> &parts, const string &mapper,
                vector<MappedRead> &mapped_reads) {
  if (mapper == "novoalign") parseMappedRead_novoAlign(parts, mapped_reads);
  else if (mapper == "bowtie") read_bowtie_output(parts, mapped_reads);
  else if (mapper == "rmap") read_rmap_output(parts, mapped_reads);
  else throw SMITHLABException("The mapper '" + mapper +\
                               "' was not recognized!");
}


/******************************************************************************
 *
 ******************************************************************************/

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

/******************************************************************************
 *                        FUNCTIONS FOR SEQUENCE I/O
 ******************************************************************************/

/***
 * \summary load sequences from either a fasta file, or a set of genomic
 *          regions and chromosome files.
 * \param target_file should reference either a fasta file or a bed file
 * \param chrom_dir   if <target_file> is a bed file, this must reference the
 *                    directory from which to load the sequences.
 * \param sequences   the sequences loaded (as strings) will be placed here
 * \param names       will be filled with the names of the sequences; either
 *                    the name as it appeared in the fasta file, or a name
 *                    with the format chr:start-end if loaded from genomic
 *                    regions and chrom. files
 * \param targets     the genomic regions associated with the sequences will be
 *                    placed here if a .bed file was provided, otherwise this
 *                    is left unchanged.
 * \param padding     pad the sequences by this much on each side. If they
 *                    came from chromosome files, the padding will be real
 *                    sequence data taken from the flanking regions; if they
 *                    came from a fasta file, the padding will be N's.
 *
 * \throw SMITHLABException if: a .bed file is given, but no valid chrom. dir.
 *                              a .fa file is given, and a chrom. dir.
 *                              the targets_file is not readable.
 * TODO this function has poor cohesion, way too specific for a library func.
 */
void
load_sequences(const string &targets_file,
               const string &chrom_dir,
               vector<string> &sequences,
               vector<string> &names,
               vector<GenomicRegion> &targets,
               const size_t padding) {

  std::ifstream in(targets_file.c_str(), std::ios::binary);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(targets_file));

  string buffer;
  getline(in, buffer);
  in.close();

  if (buffer[0] == '>') {
    if (!chrom_dir.empty()) {
      throw SMITHLABException("Provided chromosomes directory, but loading "
                              "sequences from fasta file; chromosomes "
                              "directory cannot be used.");
    }
    read_fasta_file(targets_file, names, sequences);
    if (padding != 0) {
      for (size_t i = 0; i < sequences.size(); ++i)
        sequences[i] = string(padding, 'N') + sequences[i] +\
                       string(padding, 'N');
    }
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


/******************************************************************************
 *                      FUNCTIONS FOR MAPPED-READS I/O
 ******************************************************************************/

/**
 * \brief load a set of mapped reads from a given file; accepted formats are
 *        those output by novoalign, bowtie and rmap (mapped-read format). If
 *        the format is novo-align, only uniquely-mapped reads are loaded.
 * \param reads_file
 * \param mapper        the mapper used to produce the file (i.e. the format)
 * \param mapped_reads  TODO
 * \throw SMITHLABException if the specified mapper type is not known.
 */
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
      vector<string> tokens(split_extended_whitespace_quoted(buffer));
      if ((mapper != "novoalign") ||
          ((mapper == "novoalign") && (isUniqueMapper(tokens, mapper))) ||
          ((mapper == "bowtie") && (tokens.size() > 0))) {
        parseMappedRead(tokens, mapper, mapped_reads);
      }
    }
    size_t percent_done = static_cast<size_t>(in.tellg()) * 100 / end_of_data;
    if (percent_done != current_done) {
      std::cerr << "\r" << percent_done << "% completed..." << std::flush;
      current_done = percent_done;
    }
    in.peek();
  }
  in.close();
  std::cerr << std::endl << mapped_reads.size()
            << " uniquely mapped reads are loaded!" << std::endl;
}

/**
 * \brief fill a buffer with mapped reads from an input stream; accepted
 *        formats are those output by novoalign, bowtie and rmap (mapped-read
 *        format). If the format is novo-align, only uniquely-mapped reads
 *        are loaded.
 * \param in        TODO
 * \param mapper    TODO
 * \param buffer    TODO
 */
void
fill_buffer_mapped_reads(std::ifstream &in, const string &mapper,
                         vector<MappedRead> &buffer) {
  static const size_t line_buffer_size = 10000;
  size_t i = 0;
  while (i < buffer.size() && !in.eof()) {
    char line_buffer[line_buffer_size];
    in.getline(line_buffer, line_buffer_size);
    if (in.gcount() == line_buffer_size - 1)
      throw SMITHLABException("Line too long");
    if (!is_header_line(line_buffer) && !is_track_line(line_buffer)) {
      vector<string> tokens(split_extended_whitespace_quoted(line_buffer));
      if (((mapper == "novoalign") && (isUniqueMapper(tokens, mapper))) ||
          ((mapper == "bowtie") && (tokens.size() > 0)) ||
          ((mapper == "rmap") && (tokens.size() > 0))) {
        MappedRead mr;
        parseMappedRead(tokens, mapper, mr);
        buffer[i] = mr;
        ++i;
      }
    }
  }
  if (i < buffer.size()) buffer.erase(buffer.begin() + i, buffer.end());
}

/******************************************************************************
 *                    FUNCTIONS FOR DIAGNOSTIC EVENT I/O
 ******************************************************************************/

/**
 * \brief load the diagnostic events from the given filename. The format of
 *        the file should be one sequence per line, each line is a comma
 *        separated list, where each element is the number of diagnostic
 *        events observed at that location in the given sequence.
 * \param fn            the filename to read the events from
 * \param diagEvents    the events are added to this vector. Any existing data
 *                      is cleared from the vector. diagEvents[i][j] is
 *                      fraction of diagnostic events in sequence i that fall
 *                      at location j.
 * \return the count of diagnostic events that were found
 */
double
loadDiagnosticEvents(const string &fn, vector<vector<double> > &diagEvents,
                     const float epsilon) {
  const bool DEBUG = false;
  size_t total = 0;
  static const size_t buffer_size = 100000; // TODO magic number
  ifstream in(fn.c_str());
  if (!in) throw SMITHLABException("failed to open input file " + fn);

  diagEvents.clear();
  while (!in.eof()) {
    // getting and splitting the line
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw SMITHLABException("Line too long in file: " + fn);
    vector<string> parts(smithlab::split(string(buffer), ",", false));
    if (parts.size() == 0) continue;
    if (DEBUG) {
      std::cerr << "line is: " << string(buffer) << endl;
      std::cerr << "split into " << parts.size() << " parts" << endl;
    }

    // parsing the line into a vector of counts ---------------->
    // get the raw counts, with a pseudo-count added, and remember the total
    // number of DEs in this sequence, before scaling, but after psuedo-count
    vector<double> diag_counts;
    for (size_t i = 0; i < parts.size(); ++i) {
      // convert the string to a size_t, being careful to check for vals < 0
      int t_val = atoi(parts[i].c_str());
      if (t_val < 0)
        throw SMITHLABException("line contains negative counts: " +\
                                string(buffer));
      size_t st_val = static_cast<size_t>(t_val);
      diag_counts.push_back(st_val + 1.0);
      total += st_val;
    }
    double max_peak = *std::max_element(diag_counts.begin(), diag_counts.end());

    // scale the counts up using epsilon
    double min_val = (epsilon * max_peak);
    for (size_t i = 0; i < diag_counts.size(); ++i) {
      diag_counts[i] = std::max(min_val, diag_counts[i]);
    }
    double total_di_after_scaling = std::accumulate(diag_counts.begin(),
                                                    diag_counts.end(), 0.0);
    if (DEBUG) {
      std::cerr << "after scalling:" << endl;
      for (size_t i = 0; i < diag_counts.size(); ++i) {
        std::cerr << diag_counts[i] << ",";
      }
      std::cerr << endl;
    }

    // convert the diag counts into probabilities and add them to diagEvents
    diagEvents.push_back(vector<double>());
    for (size_t i = 0; i < diag_counts.size(); ++i) {
      double frc = diag_counts[i] / static_cast<double>(total_di_after_scaling);
      diagEvents.back().push_back(frc);
    }
    if (DEBUG) {
      std::cerr << "after conversion to prob.:" << endl;
      for (size_t i = 0; i < diagEvents.back().size(); ++i)
        std::cerr << diagEvents.back()[i] << ",";
      std::cerr << endl;
    }
  }
  return total;
}


/******************************************************************************
 *                   FUNCTIONS FOR SECONDARY STRUCTURE I/O
 ******************************************************************************/

/**
 * \brief load a set of RNA secondary structures from a file. The format of the
 *        file should be one structure per line, each line/structure is a
 *        comma separated set of base-pair probabilities (floating point
 *        numbers between 0 and 1, inclusive), each one representing the
 *        probability that the base at the corresponding position within its
 *        matching sequence is paired (double stranded). The sequences
 *        themselves are not stored/represented in the file.
 * \param structure_file    filename to load from.
 * \param structures        the resultant base-pair probabilities will be
 *                          placed into this vector. Existing entries are NOT
 *                          cleared. structures[i][j] will be the base-pair
 *                          probability for the jth base in the ith sequence.
 * \throw SMITHLABException if the file cannot be read from.
 * \throw SMITHLABException if any entry v is outside of the range 0 <= v <= 1.
 */
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
      if ((d < 0) || (d > 1)) {
        stringstream ss;
        ss << "Reading structure file failed. Could not parse '" << s << "'; "
           << "resultant value was '" << d << "'. Expected floating point "
           << "value between 0 and 1. Check file format";
        throw SMITHLABException(ss.str());
      }
      record.push_back(d);
    }
    structures.push_back(record);
  }
}

/**
 * \brief write the provided secondary structures (base-pair probability
 *        vectors) to the provided file.
 * \param padding amount of padding that was added to the sequences before
 *        these structures were calculated; the values corresponding to the
 *        padding will be stripped before writing.
 * \throw SMITHLABException if file cannot be opened.
 */
void
save_structure_file(const vector<vector<double> > &sec_structure,
                    const string &outfile,
                    const size_t padding) {
  if (outfile.empty())
    throw SMITHLABException("Failed writing structure to file. Empty filename.");
  std::ofstream out(outfile.c_str());
  if (!out.good())
    throw SMITHLABException("Failed writing structure to " + outfile + "; " +\
                            "File not writable");
  save_structure_file(sec_structure, out, padding);
}


/**
 * \brief write the provided secondary structures (base-pair probability
 *        vectors) to the provided stream.
 * \param padding amount of padding that was added to the sequences before
 *        these structures were calculated; the values corresponding to the
 *        padding will be stripped before writing.
 * \throw SMITHLABException if stream is bad.
 */
void
save_structure_file(const vector<vector<double> > &sec_structure,
                    std::ostream &out,
                    const size_t padding) {
  if (!out.good())
    throw SMITHLABException("Failed writing structure to stream. Bad stream.");

  // check these before we write anything, to avoid writing a partial stream
  // if something is screwed up.
  for (size_t i = 0; i < sec_structure.size(); ++i) {
    if (sec_structure[i].size() < 2 * padding) {
      stringstream ss;
      ss << "Failed writing structure to stream, structure vector number "
         << i << " has size (" << sec_structure[i].size() << ") smaller than "
         << "area to be removed as padding (2 x " << padding << ").";
      throw SMITHLABException(ss.str());
    }
  }

  for (size_t i = 0; i < sec_structure.size(); ++i) {
    copy(sec_structure[i].begin() + padding, sec_structure[i].end() - padding,
         std::ostream_iterator<double>(out, ","));
    out << endl;
  }
}
