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
   PJU -- Jan 7 2014 -- modified get_base_pair_probability_vector to split
                        large sequences into smaller parts and fixed some
                        out-of-date documentation.

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

// integral constants
const size_t RNAUtils::SEGMENT_LENGTH;
const size_t RNAUtils::OVERLAP_LENGTH;

// initialization of non-integral constants
/// A-U(T), C-G, G-C, G-U(T), U(T)-A, U(T)-G
const int RNAUtils::RNAPair[4][4] = { { 0, 0, 0, 5 }, { 0, 0, 1, 0 }, { 0, 2, 0,
    3 }, { 6, 0, 4, 0 } };
const double RNAUtils::energy_pf = -1.0;

/**
 * \brief \todo
 * \param seq       \todo
 * \param cnstrnt   \todo
 */
double RNAUtils::get_minimum_free_energy(const string seq, const string cnstrnt) {

  vector<char> seqc(seq.size() + 1);
  copy(seq.begin(), seq.end(), seqc.begin());
  vector<char> conc(seq.size() + 1);
  copy(cnstrnt.begin(), cnstrnt.end(), conc.begin());

  MC mc;
  mc.init_pf_fold(seq.length(), energy_pf);
  return mc.getMinimumFreeEnergy(&*seqc.begin(), &*conc.begin());
}

/**
 * \brief Populate the base pair probability vector for the given set of
 *        sequences. Calculated using McCaskill's partition function.
 * \param VERBOSE  \todo
 * \param seqs     \todo
 * \param bppvs    the results are place in here. bppvs[i][j] is the
 *                 probability that the jth base in the ith sequence will
 *                 be double-stranded (paired). Any existing values in these
 *                 vectors will be cleared.
 */
void RNAUtils::get_base_pair_probability_vector(bool VERBOSE,
    const vector<string> &seqs, vector<vector<double> > &bppvs) {
  bppvs.clear();
  bppvs.resize(seqs.size());
  for (size_t i = 0; i < seqs.size(); ++i) {
    if (VERBOSE) {
      const double done = i * 100 / seqs.size();
      cerr << "\r" << "CALCULATING BASE PAIR PROBABILITIES ... (" << done
           << "% COMPLETE...)" << std::flush;
    }
    get_base_pair_probability_vector(seqs[i], bppvs[i]);
  }
  if (VERBOSE)
    cerr << "\r" << "CALCULATING BASE PAIR PROBABILITIES ... DONE"
         << endl;
}

/**
 * \brief Populate the base pair probability vector for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \param seq   the sequence to compute the bppv for.
 * \param bppv  the result is placed in here. bppv[i] is the probability that
 *              the ith base will be double-stranded (paired).
 */
void
RNAUtils::get_base_pair_probability_vector(const string seq,
                                           vector<double> &bppv) {
  string constraint(seq.size(), '.');
  get_base_pair_probability_vector(seq, constraint, bppv);
}

/**
 * \brief Populate the base pair probability vector for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \param seq     the sequence to compute the bppv for.
 * \param cnstrnt \todo
 * \param bppv    bppv[i] is the probability that the ith base will be
 *                double-stranded (paired). Any existing values in this
 *                vector will be cleared
 */
void
RNAUtils::get_base_pair_probability_vector(const string seq,
                                           const string cnstrnt,
                                           vector<double> &bppv) {

  if (seq.size() != cnstrnt.size()) {
    stringstream ss;
    ss << "Calculating base pair prob. vector failed. Reason: sequence was "
        << seq.size() << " nucleotides long, but constraint string was "
        << cnstrnt.size() << " characters";
    throw SMITHLABException(ss.str());
  }

  bppv.clear();
  bppv.resize(seq.size());

  // the PF code can't handle sequences greater than a certain size, so we
  // break seq up into overlapping pieces and do each separately. For
  // overlapping positions, we take the weighted average between the two values.
  // We can get away with this because local (~200bp) interactions dominate the
  // bpp vector.
  size_t cpos = 0;
  while(cpos < seq.length()) {
    size_t sub_len = std::min(RNAUtils::SEGMENT_LENGTH, seq.length() - cpos);

    // PF code is nasty and takes char* (non-const), and might modify them. so..
    vector<char> seqc(sub_len + 1);
    copy(seq.begin() + cpos, seq.begin() + cpos + sub_len, seqc.begin());
    vector<char> conc(sub_len + 1);
    copy(cnstrnt.begin() + cpos, cnstrnt.begin() + cpos + sub_len, conc.begin());

    vector<double> bppv_local;
    MC mc;
    mc.init_pf_fold(seq.size(), energy_pf);
    double q = mc.pf_fold(&*seqc.begin(), &*conc.begin());
    if (q > std::numeric_limits<double>::min())
      mc.getProbVector(bppv_local, sub_len);

    // incorporate this segment's bppv values into the complete answer
    assert(RNAUtils::SEGMENT_LENGTH > 2 * RNAUtils::OVERLAP_LENGTH);
    for (size_t i = 0; i < bppv_local.size(); ++i) {
      assert(cpos+i < bppv.size());
      if ((i < RNAUtils::OVERLAP_LENGTH) && (cpos != 0)) {
        double w = i / RNAUtils::OVERLAP_LENGTH;
        bppv[cpos+i] = (w * bppv_local[i]) + ((1-w) * bppv[cpos+i]);
      } else bppv[cpos+i] = bppv_local[i];
    }

    cpos = cpos + RNAUtils::SEGMENT_LENGTH - RNAUtils::OVERLAP_LENGTH;
  }
}

/**
 * \brief Populate the base pair probability matrix for the given sequence.
 *        Calculated using McCaskill's partition function.
 * \param seq  the sequence to construct the matrix for.
 * \param bppm the base-pair probability matrix; bppm[i][j] is the probability
 *        that the ith base will be paired with the jth base in <seq>.
 */
void
RNAUtils::get_base_pair_probability_matrix(const string seq,
                                           vector<vector<double> > &bppm) {
  string constraint(seq.size(), '.');
  get_base_pair_probability_matrix(seq, constraint, bppm);
}

/**
 * \brief  Populate the base pair probability matrix for the given sequence.
 *         Calculated using McCaskill's partition function.
 * \param  seq the sequence to construct the matrix for.
 * \param  cnstrnt  \todo
 * \param  bppm the base-pair probability matrix; bppm[i][j] is the probability
 *         that the ith base will be paired with the jth base in <seq>.
 *         Any existing values in these vectors will be removed by this func.
 * \return todo
 * \todo   fix this function so it works with longer sequences -- see
 *         get_base_pair_probability_vector as an example of how to do this.
 */
double
RNAUtils::get_base_pair_probability_matrix(const string seq,
                                           const string cnstrnt,
                                           vector<vector<double> > &bppm) {

  MC mc;

  if (seq.size() != cnstrnt.size()) {
    stringstream ss;
    ss << "Calculating base pair prob. matrix failed. Reason: sequence was "
        << seq.size() << " nucleotides long, but constraint string was "
        << cnstrnt.size() << " characters";
    throw SMITHLABException(ss.str());
  }

  // clean the vector
  bppm.clear();

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

