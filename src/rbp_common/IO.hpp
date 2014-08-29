/****
 Description: TODO

 --------------------

 Copyright (C) 2012
 University of Southern California,
 Emad Bahrami-Samani, Philip J. Uren, Andrew D. Smith

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

#include <iostream>
#include <vector>
#include <string>

#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

double
loadDiagnosticEvents(const std::string &fn,
                     std::vector<std::vector<double> > &diagEvents,
                     std::vector<std::vector<std::vector<double> > > &diag_values,
                     const float epsilon,
                     const double de_weight,
                     const double geo_p,
                     const size_t w);

void
load_sequences(const std::string &targets_file,
               const std::string &chrom_dir,
               std::vector<std::string> &sequences,
               std::vector<std::string> &names,
               std::vector<GenomicRegion> &targets,
               const size_t padding = 0);

void
load_mapped_reads(const std::string &reads_file,
                  const std::string &mapper,
                  std::vector<MappedRead> &mapped_reads);

void
fill_buffer_mapped_reads(std::ifstream &in, const std::string &mapper,
                         std::vector<MappedRead> &buffer);

void
load_structures(const std::string structure_file,
                std::vector<std::vector<double> > &structures);

void
save_structure_file(const std::vector<std::vector<double> > &sec_structure,
                    const std::string &outfile,
                    const size_t padding);

void
save_structure_file(const std::vector<std::vector<double> > &sec_structure,
                    std::ostream &out,
                    const size_t padding);

bool
seq_and_structure_are_consistent(const std::vector<std::string> &seqs,
                                 const std::vector<std::vector<double> > &sec_structure);


#endif
