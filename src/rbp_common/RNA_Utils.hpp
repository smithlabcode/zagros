/**
 \file RNA_Utils.hpp
 \brief Functions for computing various properties of RNA sequences and
 structures

 \authors Emad Bahrami Samani, Philip J. Uren

 \section copyright Copyright Details

 Copyright (C) 2012
 University of Southern California,
 Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

 \section license License Details

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

 \section bugs Known Bugs

 \section history Revision History

 **/

#ifndef RNA_UTILS_H
#define RNA_UTILS_H

/// scaling factor for partition function energies.

#include <vector>
#include <string>

class RNAUtils {
public:
  /*** Public static constants ***/
  /// define what bases can pair with each other
  static const int RNAPair[4][4];
  /// number of bases in the RNA alphabet
  static const size_t nbases = 4;

  static const double energy_pf;

  /*** Static functions for counting structures ***/
  static double get_minimum_free_energy(const std::string seq,
                                        const std::string cnstrnt);

  /*** static functions for determining the base pair probabilities ***/

  static void
  get_base_pair_probability_vector(bool VERBOSE,
                                   const std::vector<std::string> &seqs,
                                   std::vector<std::vector<double> > &bppvs);

  static void get_base_pair_probability_vector(const std::string seq,
                                               std::vector<double> &bppv);

  static void
  get_base_pair_probability_matrix(const std::string seq,
                                   std::vector<std::vector<double> > &bppm);

  static double get_base_pair_probability_vector(const std::string seq,
                                                 const std::string cnstrnt,
                                                 std::vector<double> &bppv);

  static double
  get_base_pair_probability_matrix(const std::string seq,
                                   const std::string constraint,
                                   std::vector<std::vector<double> > &bppm);

private:
  /// Constructor is disabled
  RNAUtils();
};

#endif
