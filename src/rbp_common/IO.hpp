/****
 Description: TODO

 --------------------

 Copyright (C) 2012
 University of Southern California,
 Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

 Authors: Emad Bahrami Samani, Philip J. Uren

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 --------------------

 Known Bugs:    None

 Revision
 History:       None
 ****/

#ifndef __INPUT_H_
#define __INPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "GenomicRegion.hpp"

struct region_less {
  bool operator()(const GenomicRegion a, const GenomicRegion b) const {
    return (a) < (b);
  }
};

class IO {
public:
  static const size_t flanking_regions_size = 20;

  static void load_sequences(std::vector<std::string> &names,
      std::vector<std::string> &sequences, std::vector<GenomicRegion> &targets,
      const std::string chrom_dir, const size_t padding,
      const std::string targets_file);

  static std::string
  print_model(const Model &model, const std::string &motif_name,
      const std::vector<GenomicRegion> &targets,
      const std::vector<std::string> &sequences,
      const std::vector<std::vector<double> > &indicators);

private:

  static void
  read_piranha_output(const std::string filename,
      std::vector<GenomicRegion> &regions);

  static void extract_regions_fasta(const std::string &dirname,
      const std::vector<GenomicRegion> &regions_in,
      std::vector<std::string> &sequences, std::vector<std::string> &names);

  static void extract_regions_chrom_fasta(const std::string &chrom_name,
      const std::string &filename, const std::vector<GenomicRegion> &regions,
      std::vector<std::string> &sequences, std::vector<std::string> &names);

  static size_t adjust_start_pos(const size_t orig_start,
      const std::string &chrom_name);

  static size_t adjust_region_size(const size_t orig_start,
      const std::string &chrom_name, const size_t orig_size);

  static void expand_regions(std::vector<GenomicRegion> &regions,
      const size_t padding);

  static void unexpand_regions(std::vector<GenomicRegion> &regions,
      const size_t padding);

  IO();
  /*** Private member constants ***/
};

#endif
