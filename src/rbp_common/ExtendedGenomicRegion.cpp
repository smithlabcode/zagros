/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2011 University of Southern California and
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

#include "ExtendedGenomicRegion.hpp"
#include "GenomicRegion.hpp"
#include "IO.hpp"

#include <cassert>
#include <fstream>
#include <sstream>
#include <tr1/unordered_map>

using std::string;
using std::vector;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::tr1::unordered_map;

#include <iostream>

unordered_map<string, chrom_id_type> ExtendedGenomicRegion::fw_table_in;
unordered_map<chrom_id_type, string> ExtendedGenomicRegion::fw_table_out;

chrom_id_type ExtendedGenomicRegion::assign_chrom(const std::string &c) {
  unordered_map<string, chrom_id_type>::const_iterator chr_id(
      fw_table_in.find(c));
  if (chr_id == fw_table_in.end()) {
    const chrom_id_type r = fw_table_in.size();
    fw_table_in[c] = r;
    fw_table_out[r] = c;
    return r;
  } else
    return chr_id->second;
}
using std::cerr;
using std::endl;

string ExtendedGenomicRegion::retrieve_chrom(chrom_id_type i) {
  unordered_map<chrom_id_type, string>::const_iterator chr_name(
      fw_table_out.find(i));
  assert(chr_name != fw_table_out.end());
  return chr_name->second;
}

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

vector<string> split_extended_whitespace_quoted(string to_split) {
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

ExtendedGenomicRegion::ExtendedGenomicRegion(string string_representation) {
  vector<string> parts(split_extended_whitespace_quoted(string_representation));

  // make sure there is the minimal required info
  if (parts.size() < 3)
    throw ExtendedGenomicRegionException(
        "Invalid string representation: " + string_representation);
  // set the chromosome name
  chrom = assign_chrom(parts[0]);

  // set the start position
  const int checkChromStart = atoi(parts[1].c_str());
  if (checkChromStart < 0)
    throw ExtendedGenomicRegionException("Invalid start: " + parts[1]);
  else
    start = static_cast<size_t>(checkChromStart);

  // set the end position
  const int checkChromEnd = atoi(parts[2].c_str());
  if (checkChromEnd < 0)
    throw ExtendedGenomicRegionException("Invalid end: " + parts[2]);
  else
    end = static_cast<size_t>(checkChromEnd);

  if (parts.size() > 3)
    name = parts[3];

  if (parts.size() > 4)
    score = atof(parts[4].c_str());

  if (parts.size() > 5)
    strand = parts[5][0];

  if (parts.size() > 6)
    sequence = parts[6];

  if (parts.size() > 7)
    extra = parts[7];

}

ExtendedGenomicRegion::ExtendedGenomicRegion(const char *s, const size_t len) {
  size_t i = 0;

  // the chrom
  while (isspace(s[i]) && i < len)
    ++i;
  size_t j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  chrom = assign_chrom(string(s + j, i - j));

  // start of the region (a positive integer)
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  start = atoi(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;

  // end of the region (a positive integer)
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  end = atoi(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;

  // name of the region
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  name = string(s + j, i - j);

  // score of the region (floating point)
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  end = atoi(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;
  score = atof(s + j);

  // strand
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  end = atoi(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;
  strand = *(s + j);

  //sequence
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  sequence = string(s + j, i - j);

  //extra
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  extra = string(s + j, i - j);
}

string ExtendedGenomicRegion::tostring() const {
  std::ostringstream s;
  s << retrieve_chrom(chrom) << "\t" << start << "\t" << end;
  if (!name.empty())
    s << "\t" << name << "\t" << score << "\t" << strand;
  if (!sequence.empty())
    s << "\t" << sequence;
  if (!extra.empty())
    s << "\t" << extra;
  return s.str();
}

std::ostream&
operator<<(std::ostream& s, const ExtendedGenomicRegion& region) {
  return s << region.tostring();
}

std::istream&
operator>>(std::istream& s, ExtendedGenomicRegion& region) {
  string chrom, name, sequence, extra;
  size_t start = 0ul, end = 0ul;
  double score = 0.0;
  char strand = '\0';

  if (s >> chrom >> start >> end >> name >> score >> strand >> sequence
      >> extra)
    region = ExtendedGenomicRegion(
        chrom, start, end, name, score, strand, sequence, extra);
  else
    region = ExtendedGenomicRegion();

  char c;
  while ((c = s.get()) != '\n' && s)
    ;

  if (c != '\n')
    s.setstate(std::ios::badbit);

  if (s.eof())
    s.setstate(std::ios::badbit);

  return s;
}

bool ExtendedGenomicRegion::contains(const ExtendedGenomicRegion& other) const {
  return chrom == other.chrom && start <= other.start && other.end <= end;
}

bool ExtendedGenomicRegion::overlaps(const ExtendedGenomicRegion& other) const {
  return chrom == other.chrom
      && ((start < other.end && other.end <= end)
          || (start <= other.start && other.start < end)
          || other.contains(*this));
}

size_t ExtendedGenomicRegion::distance(
    const ExtendedGenomicRegion& other) const {
  if (chrom != other.chrom)
    return std::numeric_limits<size_t>::max();
  else if (overlaps(other) || other.overlaps(*this))
    return 0;
  else
    return (end < other.start) ? other.start - end + 1 : start - other.end + 1;
}

bool ExtendedGenomicRegion::operator<(const ExtendedGenomicRegion& rhs) const {
  return ((chrom == rhs.chrom
      && (end < rhs.end
          || (end == rhs.end
              && (start < rhs.start
                  || (start == rhs.start && (strand < rhs.strand
                  // || (strand == rhs.strand && name < rhs.name)
                      )))))) || get_chrom() < rhs.get_chrom());
}

bool ExtendedGenomicRegion::less1(const ExtendedGenomicRegion& rhs) const {
  return ((chrom == rhs.chrom
      && (start < rhs.start
          || (start == rhs.start
              && (end < rhs.end || (end == rhs.end && (strand < rhs.strand
              // || (strand == rhs.strand && name < rhs.name)
                  )))))) || get_chrom() < rhs.get_chrom());
}

bool ExtendedGenomicRegion::operator<=(const ExtendedGenomicRegion& rhs) const {
  return !(rhs < *this);
}

bool ExtendedGenomicRegion::operator==(const ExtendedGenomicRegion& rhs) const {
  return (chrom == rhs.chrom && start == rhs.start && end == rhs.end
      && name == rhs.name && score == rhs.score && strand == rhs.strand);
}

bool ExtendedGenomicRegion::operator!=(const ExtendedGenomicRegion& rhs) const {
  return (chrom != rhs.chrom || start != rhs.start || end != rhs.end
      || name != rhs.name || score != rhs.score || strand != rhs.strand);
}

void separate_extended_chromosomes(const vector<ExtendedGenomicRegion>& regions,
    vector<vector<ExtendedGenomicRegion> >& separated_by_chrom) {
  typedef unordered_map<chrom_id_type, vector<ExtendedGenomicRegion> > Separator;
  Separator separator;
  for (vector<ExtendedGenomicRegion>::const_iterator i = regions.begin();
      i != regions.end(); ++i) {
    const chrom_id_type the_chrom(i->chrom);
    if (separator.find(the_chrom) == separator.end())
      separator[the_chrom] = vector<ExtendedGenomicRegion>();
    separator[the_chrom].push_back(*i);
  }
  separated_by_chrom.clear();
  for (Separator::iterator i = separator.begin(); i != separator.end(); ++i)
    separated_by_chrom.push_back(i->second);
}

void ReadExtendedBEDFile(string filename,
    vector<ExtendedGenomicRegion> &the_regions) {
  static const size_t buffer_size = 10000; // Magic

  // open and check the file
  std::ifstream in(filename.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + filename);
  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    the_regions.push_back(ExtendedGenomicRegion(buffer));
    in.peek();
  }
  in.close();
}


