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
#define energy_pf -1

#include <vector>
#include <string>

class RNAUtils {
public:
  /*** Public static constants ***/
  /// define what bases can pair with each other
  static const int RNAPair[4][4];
  /// number of bases in the RNA alphabet
  static const size_t nbases = 4;
  static double expLoopTable[31][31][7][7][4][4][4][4];

  static inline size_t base2int_RNA(char c) {
    switch (c) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    case 'U':
      return 3;
    case 'a':
      return 0;
    case 'c':
      return 1;
    case 'g':
      return 2;
    case 't':
      return 3;
    case 'u':
      return 3;
    default:
      return 4;
    }
  }

  /*** Static functions for counting structures ***/
  static double getMinimumFreeEnergy(std::string seq, std::string cnstrnt);
  static std::vector<std::vector<std::vector<double> > >
  calculateWeightedNumberOfStructures(const std::vector<std::string> &seqs,
      const std::vector<std::vector<std::vector<double> > > &Ps);
  static std::vector<std::vector<double> >
  calculateWeightedNumberOfStructures(const std::string &sequence,
      const std::vector<std::vector<double> > &P);

  static std::vector<std::vector<double> >
  calculateNumberOfStructures(const std::string &sequence);
  static std::vector<std::vector<std::vector<double> > >
  calculateNumberOfStructures(const std::vector<std::string> &seqs);

  /*** static functions for counting the number of (un)paired bases ***/
  static std::vector<double> calculateUnpairedState(const std::string &seq,
      const std::vector<double> &basePairProbs, double weighteNumberSecStr);
  static std::vector<double> calculatePairedState(const std::string &seq,
      const std::vector<double> &basePairProbs, double weighteNumberSecStr);
  static std::vector<std::vector<double> >
  calculateUnpairedState(const std::vector<std::string> &seqs,
      const std::vector<std::vector<double> > &basePairProbs,
      const std::vector<double> &wnSecStr);
  static std::vector<std::vector<double> >
  calculatePairedState(const std::vector<std::string> &seqs,
      const std::vector<std::vector<double> > &basePairProbs,
      const std::vector<double> &wnSecStr);

  /*** static functions for counting bases in sequences ***/
  static std::vector<std::vector<size_t> >
  baseOccurrenceIndicatorMatrix(std::string sequence);
  static std::vector<std::vector<std::vector<size_t> > >
  baseOccurrenceIndicatorMatrix(std::vector<std::string> sequences);

  /*** static functions for determining the base pair probabilities ***/
  static std::vector<double> getBasePairProbabilityVector(std::string seq);
  static std::vector<std::vector<double> >
  getBasePairProbabilityMatrix(std::string seq);
  static std::vector<double> getBasePairProbabilityVector(std::string seq,
      std::string cnstrnt);
  static std::vector<std::vector<double> >
  getBasePairProbabilityMatrix(std::string seq, std::string constraint);

  static void fillExpLoopTable();

  static std::string getReverse(std::string sequence);

  static std::string getBaseComplement(std::string base);

  static std::string getComplement(std::string sequence);

  static std::string getReverseComplement(std::string sequence);

private:
  /// Constructor is disabled
  RNAUtils();
};

#endif
