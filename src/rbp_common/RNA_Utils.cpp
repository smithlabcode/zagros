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
#include <sstream>
#include <fstream>
#include <bitset>
#include <limits>
#include <float.h> 

#include "RNA_Utils.hpp"
#include "smithlab_utils.hpp"
#include "part_func.hpp"

using std::string;
using std::vector;
using std::stringstream;
using std::ifstream;
using std::ostream;
using std::ostream_iterator;
using std::ofstream;
using std::cerr;
using std::endl;

/// A-U(T), C-G, G-C, G-U(T), U(T)-A, U(T)-G
const int RNAUtils::RNAPair[4][4] = { { 0, 0, 0, 5 }, { 0, 0, 1, 0 }, { 0, 2, 0,
    3 }, { 6, 0, 4, 0 } };

double RNAUtils::expLoopTable[31][31][7][7][4][4][4][4];

string RNAUtils::getReverse(string sequence) {
  string s = "";
  for (size_t i = 0; i < sequence.length(); ++i)
    s = sequence.substr(i, 1) + s;
  return s;
}

string RNAUtils::getBaseComplement(string base) {
  if (base == "A")
    return "T";
  else if (base == "C")
    return "G";
  else if (base == "G")
    return "C";
  else
    return "A";
}

string RNAUtils::getComplement(string sequence) {
  string s = "";
  for (size_t i = 0; i < sequence.length(); ++i)
    s = s + getBaseComplement(sequence.substr(i, 1));
  return s;
}

string RNAUtils::getReverseComplement(string sequence) {
  string s = "";
  for (size_t i = 0; i < sequence.length(); ++i)
    s = getBaseComplement(sequence.substr(i, 1)) + s;
  return s;
}

void RNAUtils::fillExpLoopTable() {

  MC mc;
  mc.init_pf_fold(120, energy_pf);
  for (size_t i = 0; i < 31; i++) {
    for (size_t j = 0; j < 31; j++) {
      for (size_t k = 0; k < 7; k++) {
        for (size_t l = 0; l < 7; l++) {
          for (size_t a1 = 0; a1 < 4; a1++) {
            for (size_t a2 = 0; a2 < 4; a2++) {
              for (size_t a3 = 0; a3 < 4; a3++) {
                for (size_t a4 = 0; a4 < 4; a4++) {
                  expLoopTable[i][j][k][l][a1][a2][a3][a4] = mc.expLoopEnergy(
                      i, j, k, l, a1, a2, a3, a4);
                }
              }
            }
          }
        }
      }
    }
  }
}

double RNAUtils::getMinimumFreeEnergy(string seq, string cnstrnt) {

  std::vector<char> seqc(seq.size() + 1);
  std::copy(seq.begin(), seq.end(), seqc.begin());
  std::vector<char> conc(seq.size() + 1);
  std::copy(cnstrnt.begin(), cnstrnt.end(), conc.begin());

  MC mc;
  mc.init_pf_fold(seq.length(), energy_pf);
  return mc.getMinimumFreeEnergy(&*seqc.begin(), &*conc.begin());
}

/**
 * \brief calculate the number of times base i appears as single stranded in
 *        any structure formed from the given sequence
 * \param seq the input sequence to use
 * \param basePairProbs vector where ith entry gives the probability that
 *        the ith base is paired
 * \return alpha where alpha[i] is the normalised count of the number of times
 *         that base i appears unpaired in any structure.
 */
vector<double> RNAUtils::calculateUnpairedState(const string &seq,
    const vector<double> &basePairProbs, double weighteNumberSecStr) {
  const size_t L = seq.size();
  vector<double> alpha(L, 0);

  for (size_t i = 0; i < L; i++) {
    alpha[i] = (1 - basePairProbs[i]) * weighteNumberSecStr;
  }

  double max = 0;
  for (size_t i = 0; i < L; i++)
    if (alpha[i] > max)
      max = alpha[i];
  for (size_t i = 0; i < L; i++)
    alpha[i] = alpha[i] / (max);

  return alpha;
}

/**
 * \brief calculate the number of times base i appears as double stranded in
 *        any structure formed from the given sequence
 * \param seq the input sequence to use
 * \param basePairProbs vector where ith entry is the probability that ith
 *                      nucleotide is paired
 * \return alpha where alpha[i] is the normalised count of the number of times
 *         that base i appears paired in any structure.
 */
vector<double> RNAUtils::calculatePairedState(const string &seq,
    const vector<double> &basePairProbs, double weighteNumberSecStr) {
  const size_t L = seq.size();
  vector<double> gamma(L, 0);
  for (size_t i = 0; i < L; i++) {
    gamma[i] = basePairProbs[i] * weighteNumberSecStr;
  }

  double max = 0;
  for (size_t i = 0; i < L; i++)
    if (gamma[i] > max)
      max = gamma[i];
  for (size_t i = 0; i < L; i++)
    gamma[i] = gamma[i] / (max);

  return gamma;
}

/**
 * \brief calculate the normalised count of times each base appears
 *        single-stranded in any structure
 * \param seqs vector of sequences
 * \return matrix M where m[i][j] is the normalised count of the number of
 *         times that base j appears unpaired in any structure from sequence i
 */
vector<vector<double> > RNAUtils::calculateUnpairedState(
    const vector<string> &seqs, const vector<vector<double> > &basePairProbs,
    const vector<double> &wnSecStr) {
  vector<vector<double> > res;
  for (size_t i = 0; i < seqs.size(); i++) {
    res.push_back(
        calculateUnpairedState(seqs[i], basePairProbs[i], wnSecStr[i]));
  }
  return res;
}

/**
 * \brief calculate the normalised count of times each base appears
 *        double-stranded in any structure
 * \param seqs vector of sequences
 * \return matrix M where m[i][j] is the normalised count of the number of
 *         times that the jth base appears paired in any structure from
 *         sequence i
 */
vector<vector<double> > RNAUtils::calculatePairedState(
    const vector<string> &seqs, const vector<vector<double> > &basePairProbs,
    const vector<double> &wnSecStr) {
  vector<vector<double> > res;
  for (size_t i = 0; i < seqs.size(); i++) {
    res.push_back(calculatePairedState(seqs[i], basePairProbs[i], wnSecStr[i]));
  }
  return res;
}

/**
 * \brief TODO
 */
vector<vector<vector<double> > > RNAUtils::calculateWeightedNumberOfStructures(
    const vector<string> &seqs, const vector<vector<vector<double> > > &Ps) {
  vector<vector<vector<double> > > res;
  for (size_t i = 0; i < seqs.size(); i++) {
    res.push_back(calculateWeightedNumberOfStructures(seqs[i], Ps[i]));
  }
  return res;
}

/**
 * \brief For a given sequence, calculate a weighted count of the number of
 *        structures possible
 * \param P The base pair probability matrix P[i][j] = prob of ith base in the
 *          sequence being paired with the jth base
 * \return An L x L matrix M where M[i][j] is the number of structures that
 *         can be produced from the subsequence from i to j, weighted by
 *         the probability of each structure occuring
 */
vector<vector<double> > RNAUtils::calculateWeightedNumberOfStructures(
    const string &sequence, const vector<vector<double> > &P) {
  const size_t L = sequence.size();
  vector<vector<double> > res;
  res.resize(L);

  // make the LxL matrix and fill with zeros
  for (size_t i = 0; i < L; i++) {
    res[i].resize(L);
    for (size_t j = 0; j < L; j++) {
      res[i][j] = 0;
    }
  }

  // set the first four diagonals through the matrix to be 1.
  // i.e. only 1 structure is possible here because we don't allow internal
  // loops of less than 4 nucleotides.
  for (size_t l = 1; l < 5; l++) {
    for (size_t i = 0; i <= L - l; i++) {
      size_t j = i + l - 1;
      res[i][j] = 1;
    }
  }

  for (size_t l = 5; l <= L; l++) {
    for (size_t i = 0; i <= L - l; i++) {
      size_t j = i + l - 1;
      string sseq = sequence.substr(i, l);
      string cnstrnt = "";
      cnstrnt.resize(l, '.');
      vector<double> T = getBasePairProbabilityVector(sseq, cnstrnt);

      res[i][j] += res[i][j - 1] * (1 - T[l - 1]);
      if (res[i][j] < 1)
        res[i][j] = 1;

      for (size_t h = i + 1; h < j - 4; h++) {
        if (RNAPair[base2int_RNA(sequence[h])][base2int_RNA(sequence[j])] != 0)
          res[i][j] += res[i][h - 1] * res[h + 1][j - 1] * P[h][j];
      }
      if (RNAPair[base2int_RNA(sequence[i])][base2int_RNA(sequence[j])] != 0)
        res[i][j] += res[i + 1][j - 1] * P[i][j];
    }
  }
  return res;
}

/**
 * \brief TODO
 */
std::vector<std::vector<double> > RNAUtils::calculateNumberOfStructures(
    const std::string &sequence) {
  vector<vector<double> > P;
  const size_t L = sequence.size();
  for (size_t l = 0; l < L; l++) {
    P.push_back(vector<double>(L, 1));
  }
  return calculateWeightedNumberOfStructures(sequence, P);
}

/**
 * \brief TODO
 */
std::vector<std::vector<std::vector<double> > > RNAUtils::calculateNumberOfStructures(
    const std::vector<std::string> &seqs) {
  vector<vector<vector<double> > > res;
  for (size_t i = 0; i < seqs.size(); i++) {
    res.push_back(calculateNumberOfStructures(seqs[i]));
  }
  return res;
}

/**
 * \brief Calculate a matrix of indicator variables specifying which bases
 *        appear at which locations in the given sequence
 * \return matrix M where M[b][l] = 1 if the sequence contains base b
 *         at location l
 */
vector<vector<size_t> > RNAUtils::baseOccurrenceIndicatorMatrix(
    std::string sequence) {
  const size_t L = sequence.size();

  vector<vector<size_t> > nu(nbases, vector<size_t>(L, 0));

  for (size_t j = 0; j < sequence.size(); j++) {
    for (size_t b = 0; b < nbases; b++) {
      nu[b][j] = 0;
    }
    nu[RNAUtils::base2int_RNA(sequence[j])][j] = 1;
  }

  return nu;
}

/**
 * \brief Calculate a matrix of indicator variables specifying which bases
 *        appear at which locations in the given sequences
 * \return matrix M where M[b][l][n] = 1 if the nth sequence contains base b
 *         at location l
 */
vector<vector<vector<size_t> > > RNAUtils::baseOccurrenceIndicatorMatrix(
    std::vector<std::string> sequences) {
  const size_t N = sequences.size();
  size_t L = 0;
  for (size_t i = 0; i < sequences.size(); i++)
    if (sequences[i].size() > L)
      L = sequences[i].size();
  vector<vector<vector<size_t> > > nu(
      nbases, vector<vector<size_t> >(L, vector<size_t>(N, 0)));
  for (size_t n = 0; n < N; n++) {
    string s1 = sequences[n];
    for (size_t j = 0; j < s1.size(); j++) {
      for (size_t b = 0; b < nbases; b++) {
        nu[b][j][n] = 0;
      }
      nu[RNAUtils::base2int_RNA(s1[j])][j][n] = 1;
    }
  }
  return nu;
}

/**
 * \brief Return the base pair probability vector for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return vector V where V[i] is the probability that the ith base will
 *         be double-stranded (paired)
 */
std::vector<double> RNAUtils::getBasePairProbabilityVector(std::string seq) {
  string constraint(seq.size(), '.');
  return getBasePairProbabilityVector(seq, constraint);
}

/**
 * \brief Return the base pair probability vector for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return vector V where V[i] is the probability that the ith base will
 *         be double-stranded (paired)
 * \todo cache results
 * \todo make MC object static member
 */
vector<double> RNAUtils::getBasePairProbabilityVector(std::string seq,
    std::string cnstrnt) {
  const size_t L = seq.size();

  if (L != cnstrnt.size()) {
    stringstream ss;
    ss << "Calculating base pair prob. vector failed. Reason: sequence was "
        << L << " nucleotides long, but constraint string was "
        << cnstrnt.size() << " characters";
    throw SMITHLABException(ss.str());
  }

  // PF code is nasty and takes char* (non-const), and might modify them.. so..
  std::vector<char> seqc(seq.size() + 1);
  std::copy(seq.begin(), seq.end(), seqc.begin());
  std::vector<char> conc(seq.size() + 1);
  std::copy(cnstrnt.begin(), cnstrnt.end(), conc.begin());

  vector<double> basePairProbs(L, 0);
  MC mc;
  mc.init_pf_fold(L, energy_pf);
  float q = mc.pf_fold(&*seqc.begin(), &*conc.begin());
  if (q <= FLT_MIN)
    return basePairProbs;
  else {
    mc.getProbVector(basePairProbs);
    return basePairProbs;
  }
}

/**
 * \brief Return the base pair probability matrix for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return matrix M where M[i][j] is the probability that the ith base will
 *         be paired with the jth base.
 */
std::vector<std::vector<double> > RNAUtils::getBasePairProbabilityMatrix(
    std::string seq) {
  string constraint(seq.size(), '.');
  return getBasePairProbabilityMatrix(seq, constraint);
}

/**
 * \brief Return the base pair probability matrix for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \return matrix M where M[i][j] is the probability that the ith base will
 *         be paired with the jth base.
 * \todo cache results
 * \todo make MC object static member
 */
std::vector<std::vector<double> > RNAUtils::getBasePairProbabilityMatrix(
    std::string seq, std::string cnstrnt) {
  MC mc;
  const size_t L = seq.size();

  if (L != cnstrnt.size()) {
    stringstream ss;
    ss << "Calculating base pair prob. matrix failed. Reason: sequence was "
        << L << " nucleotides long, but constraint string was "
        << cnstrnt.size() << " characters";
    throw SMITHLABException(ss.str());
  }

  // PF code is nasty and takes char* (non-const), and might modify them.. so..
  std::vector<char> seqc(seq.size() + 1);
  std::copy(seq.begin(), seq.end(), seqc.begin());
  std::vector<char> conc(seq.size() + 1);
  std::copy(cnstrnt.begin(), cnstrnt.end(), conc.begin());

  mc.init_pf_fold(L, energy_pf);
  vector<vector<double> > P(L, vector<double>(L, 0));

  float q = mc.pf_fold(&*seqc.begin(), &*conc.begin());
  if (q <= FLT_MIN)
    return P;
  else {
    mc.getProbMatrix(P, L);
    return P;
  }
}

