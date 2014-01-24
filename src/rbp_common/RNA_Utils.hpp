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

#include <vector>
#include <string>

#include "RNG.hpp"


// TODO -- PJU: probably no need for this to be a class; a namespace would
//              be better since everything is public and static and the c'tor
//              is disabled...
class RNAUtils {
public:
  /*** Public static constants ***/
  /// define what bases can pair with each other
  static const int RNAPair[4][4];

  /// number of bases in the RNA alphabet
  static const size_t RNA_ALPHABET_SIZE = 4;

  /// scaling factor for partition function energies.
  static const double energy_pf;

  // sequences that are longer than this will be broken into overlapping
  // tiles for determining their base-pair probabilities.
  static const size_t SEGMENT_LENGTH = 250;

  // large sequences that are broken into tiles will have this much overlap
  // between each tile.
  static const size_t OVERLAP_LENGTH = 50;

  // default background distribution for nucleotide occurrences in RNA
  const static double DEFAULT_BACKGROUND[RNA_ALPHABET_SIZE];

  /*** indexing ***/
  static inline size_t
  base2int(char c) {
    switch(c) {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    case 'U' : return 3;
    case 'a' : return 0;
    case 'c' : return 1;
    case 'g' : return 2;
    case 't' : return 3;
    case 'u' : return 3;
    default  : return 4;
    }
  }

  /*** Static functions for generating random RNA sequences ***/
  static char sampleNuc(const Runif &rng);
  static char sampleNucFromNucDist(const std::vector<double> &dist,
                                   const Runif &rng);
  static std::string sampleSeq(const size_t length, const Runif &rng);
  static std::string sampleSeqFromNucDist(const std::vector<double> &dist,
                                          const size_t length,
                                          const Runif &rng);
  static std::string sampleSequenceFromPWM(const std::vector<std::vector<double> > pwm,
                                           const Runif &rng);

  /*** Static functions for counting structures ***/
  static double get_minimum_free_energy(const std::string seq,
                                        const std::string cnstrnt);

  /*** static functions for determining the base pair probabilities ***/
  static void
  get_base_pair_probability_vector(bool VERBOSE,
                                   const std::vector<std::string> &seqs,
                                   std::vector<std::vector<double> > &bppvs);

  static void
  get_base_pair_probability_vector(const std::string seq,
                                   std::vector<double> &bppv);

  static void
  get_base_pair_probability_matrix(const std::string seq,
                                   std::vector<std::vector<double> > &bppm);

  static void
  get_base_pair_probability_vector(const std::string seq,
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
