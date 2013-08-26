/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
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

#ifndef EXTENDED_GENOMIC_REGION_HPP
#define EXTENDED_GENOMIC_REGION_HPP

#include "smithlab_utils.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <tr1/unordered_map>
#include <limits>

typedef unsigned chrom_id_type;

class ExtendedGenomicRegion {
public:
  ExtendedGenomicRegion() :
      chrom(assign_chrom("(NULL)")), name("X"), start(0), end(0), score(0), strand(
          '+'), sequence("Y"), extra("Z") {
  }
  void swap(ExtendedGenomicRegion &rhs) {
    std::swap(chrom, rhs.chrom);
    std::swap(name, rhs.name);
    std::swap(start, rhs.start);
    std::swap(end, rhs.end);
    std::swap(score, rhs.score);
    std::swap(strand, rhs.strand);
    std::swap(sequence, rhs.sequence);
    std::swap(extra, rhs.extra);
  }
  ExtendedGenomicRegion(const ExtendedGenomicRegion &other) :
      chrom(other.chrom), name(other.name), start(other.start), end(other.end), score(
          other.score), strand(other.strand), sequence(other.sequence), extra(
          other.extra) {
  }
  ExtendedGenomicRegion& operator=(const ExtendedGenomicRegion& rhs) {
    ExtendedGenomicRegion tmp(rhs);
    swap(tmp);
    return *this;
  }

  // Other constructors
  ExtendedGenomicRegion(std::string c, size_t sta, size_t e, std::string n,
      float sc, char str, std::string s, std::string ext) :
      chrom(assign_chrom(c)), name(n), start(sta), end(e), score(sc), strand(
          str), sequence(s), extra(ext) {
  }
  ExtendedGenomicRegion(std::string c, size_t sta, size_t e) :
      chrom(assign_chrom(c)), name(std::string("X")), start(sta), end(e), score(
          0.0), strand('+'), sequence("Y"), extra("Z") {
  }
  explicit ExtendedGenomicRegion(std::string string_representation);
  ExtendedGenomicRegion(const char *s, const size_t len);
  std::string tostring() const;

  // accessors
  std::string get_chrom() const {
    return retrieve_chrom(chrom);
  }
  size_t get_start() const {
    return start;
  }
  size_t get_end() const {
    return end;
  }
  size_t get_width() const {
    return (end > start) ? end - start : 0;
  }
  std::string get_name() const {
    return name;
  }
  float get_score() const {
    return score;
  }
  char get_strand() const {
    return strand;
  }
  bool pos_strand() const {
    return (strand == '+');
  }
  bool neg_strand() const {
    return (strand == '-');
  }
  std::string get_sequence() const {
    return sequence;
  }
  std::string get_extra() const {
    return extra;
  }

  // mutators
  void set_chrom(const std::string& new_chrom) {
    chrom = assign_chrom(new_chrom);
  }
  void set_start(size_t new_start) {
    start = new_start;
  }
  void set_end(size_t new_end) {
    end = new_end;
  }
  void set_name(const std::string& n) {
    name = n;
  }
  void set_score(float s) {
    score = s;
  }
  void set_strand(char s) {
    strand = s;
  }
  void set_sequence(const std::string& s) {
    sequence = s;
  }
  void set_extra(const std::string& e) {
    extra = e;
  }

  // comparison functions
  bool contains(const ExtendedGenomicRegion& other) const;
  bool overlaps(const ExtendedGenomicRegion& other) const;
  size_t distance(const ExtendedGenomicRegion& other) const;
  bool operator<(const ExtendedGenomicRegion& rhs) const;
  bool less1(const ExtendedGenomicRegion& rhs) const;
  bool operator<=(const ExtendedGenomicRegion& rhs) const;
  bool operator!=(const ExtendedGenomicRegion& rhs) const;
  bool operator==(const ExtendedGenomicRegion& rhs) const;

  bool same_chrom(const ExtendedGenomicRegion &other) const {
    return chrom == other.chrom;
  }

  friend void
  separate_extended_chromosomes(
      const std::vector<ExtendedGenomicRegion>& regions,
      std::vector<std::vector<ExtendedGenomicRegion> >& separated_by_chrom);

private:

  static chrom_id_type assign_chrom(const std::string &c);
  static std::string retrieve_chrom(chrom_id_type i);

  static std::tr1::unordered_map<std::string, chrom_id_type> fw_table_in;
  static std::tr1::unordered_map<chrom_id_type, std::string> fw_table_out;
  // std::string chrom;
  chrom_id_type chrom;
  std::string name;
  size_t start;
  size_t end;
  float score;
  char strand;
  std::string sequence;
  std::string extra;
};

class ExtendedGenomicRegionException: public SMITHLABException {
public:
  ExtendedGenomicRegionException(std::string s = std::string()) :
      SMITHLABException(s) {
  }
};

size_t
convertString(const std::string& s);

std::string
convertSizet(const size_t number);

std::vector<std::string>
split_extended_whitespace_quoted(std::string to_split);

void
ReadExtendedBEDFile(std::string filename,
    std::vector<ExtendedGenomicRegion> &regions);

void
read_rmap_output(std::string filename,
    std::vector<ExtendedGenomicRegion> &regions);

void
read_novoalign_output(std::string filename,
    std::vector<ExtendedGenomicRegion> &regions);

void
read_piranha_output(std::string filename,
    std::vector<ExtendedGenomicRegion> &regions);

void
read_bowtie_output(std::string filename,
    std::vector<ExtendedGenomicRegion> &regions);

template<class T> void WriteExtendedBEDFile(const std::string filename,
    const std::vector<std::vector<T> > &regions, std::string track_name = "") {
  std::ofstream out(filename.c_str());
  if (track_name.length() > 0)
    out << "track name=" << track_name << std::endl;
  for (typename std::vector<std::vector<T> >::const_iterator i =
      regions.begin(); i != regions.end(); ++i)
    std::copy(i->begin(), i->end(), std::ostream_iterator<T>(out, "\n"));
  out.close();
}

template<class T> void WriteExtendedBEDFile(const std::string filename,
    const std::vector<T> &regions, std::string track_name = "") {
  std::ofstream out(filename.c_str());
  if (track_name.length() > 0)
    out << "track name=" << track_name << std::endl;
  std::copy(
      regions.begin(), regions.end(), std::ostream_iterator<T>(out, "\n"));
  out.close();
}

#endif
