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

#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include "RNA_Utils.hpp"
#include "smithlab_utils.hpp"
#include "Part_Func.hpp"

using std::string;
using std::vector;
using std::stringstream;
using std::cerr;
using std::endl;

/// A-U(T), C-G, G-C, G-U(T), U(T)-A, U(T)-G
const int RNAUtils::RNAPair[4][4] = { { 0, 0, 0, 5 }, { 0, 0, 1, 0 }, { 0, 2, 0,
    3 }, { 6, 0, 4, 0 } };

double RNAUtils::getMinimumFreeEnergy(const string seq, const string cnstrnt) {

  vector<char> seqc(seq.size() + 1);
  copy(seq.begin(), seq.end(), seqc.begin());
  vector<char> conc(seq.size() + 1);
  copy(cnstrnt.begin(), cnstrnt.end(), conc.begin());

  MC mc;
  mc.init_pf_fold(seq.length(), energy_pf);
  return mc.getMinimumFreeEnergy(&*seqc.begin(), &*conc.begin());
}

/**
 * \brief Return the base pair probability vector for the given set of sequences
 *        Calculated using McCaskill's partition function.
 * \return matrix V where V[i][j] is the probability that the jth base in
 *         ith sequence will be double-stranded (paired)
 */
void RNAUtils::getBasePairProbabilityVector(bool VERBOSE,
    const vector<string> &seqs, vector<vector<double> > &bppvs) {
  bppvs.clear();
  bppvs.resize(seqs.size());
  for (size_t i = 0; i < seqs.size(); ++i) {
    if (VERBOSE)
      cerr << "\r" << i * 100 / seqs.size() << "% completed..."
          << std::flush;
    getBasePairProbabilityVector(seqs[i], bppvs[i]);
  }
  if (VERBOSE)
    cerr << "\r100% completed..."
        << endl;
}

/**
 * \brief Return the base pair probability vector for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return vector V where V[i] is the probability that the ith base will
 *         be double-stranded (paired)
 */
void RNAUtils::getBasePairProbabilityVector(const string seq,
    vector<double> &bppv) {
  string constraint(seq.size(), '.');
  getBasePairProbabilityVector(seq, constraint, bppv);
}

/**
 * \brief Return the base pair probability vector for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return vector V where V[i] is the probability that the ith base will
 *         be double-stranded (paired)
 * \todo cache results
 * \todo make MC object static member
 */
double RNAUtils::getBasePairProbabilityVector(const string seq,
    const string cnstrnt, vector<double> &bppv) {

  if (seq.size() != cnstrnt.size()) {
    stringstream ss;
    ss << "Calculating base pair prob. vector failed. Reason: sequence was "
        << seq.size() << " nucleotides long, but constraint string was "
        << cnstrnt.size() << " characters";
    throw SMITHLABException(ss.str());
  }

  // PF code is nasty and takes char* (non-const), and might modify them.. so..
  vector<char> seqc(seq.size() + 1);
  copy(seq.begin(), seq.end(), seqc.begin());
  vector<char> conc(seq.size() + 1);
  copy(cnstrnt.begin(), cnstrnt.end(), conc.begin());

  MC mc;
  mc.init_pf_fold(seq.size(), energy_pf);
  double q = mc.pf_fold(&*seqc.begin(), &*conc.begin());
  if (q > std::numeric_limits<double>::min())
    mc.getProbVector(bppv, seq.length());
  return q;
}

/**
 * \brief Return the base pair probability matrix for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return matrix M where M[i][j] is the probability that the ith base will
 *         be paired with the jth base.
 */
void RNAUtils::getBasePairProbabilityMatrix(const string seq,
    vector<vector<double> > &bppm) {
  string constraint(seq.size(), '.');
  getBasePairProbabilityMatrix(seq, constraint, bppm);
}

/**
 * \brief Return the base pair probability matrix for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return matrix M where M[i][j] is the probability that the ith base will
 *         be paired with the jth base.
 * \todo cache results
 * \todo make MC object static member
 */
double RNAUtils::getBasePairProbabilityMatrix(const string seq,
    const string cnstrnt, vector<vector<double> > &bppm) {

  MC mc;

  if (seq.size() != cnstrnt.size()) {
    stringstream ss;
    ss << "Calculating base pair prob. matrix failed. Reason: sequence was "
        << seq.size() << " nucleotides long, but constraint string was "
        << cnstrnt.size() << " characters";
    throw SMITHLABException(ss.str());
  }

  // PF code is nasty and takes char* (non-const), and might modify them.. so..
  vector<char> seqc(seq.size() + 1);
  copy(seq.begin(), seq.end(), seqc.begin());
  vector<char> conc(seq.size() + 1);
  copy(cnstrnt.begin(), cnstrnt.end(), conc.begin());

  mc.init_pf_fold(seq.size(), energy_pf);

  double q = mc.pf_fold(&*seqc.begin(), &*conc.begin());
  if (q > std::numeric_limits<double>::min())
    mc.getProbMatrix(bppm, seq.length());
  return q;
}

