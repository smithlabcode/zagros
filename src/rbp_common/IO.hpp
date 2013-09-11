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

#include <iostream>
#include <vector>
#include <string>

#include "GenomicRegion.hpp"


void 
load_sequences(const std::string &chrom_dir, 
	       const size_t padding,
	       const std::string &targets_file,
	       std::vector<std::string> &names,
	       std::vector<std::string> &sequences, 
	       std::vector<GenomicRegion> &targets);


#endif
