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
const double RNAUtils::DEFAULT_BACKGROUND[RNAUtils::RNA_ALPHABET_SIZE] =\
    {0.3,  0.2, 0.2,  0.3};
const double RNAUtils::energy_pf = -1.0;

/**
 * \brief sample a single nucleotide from a distribution
 * \param dist  the distribution to use; dist[base2int('A')] gives the
 *              probability that an 'A' will be sampled, and so on.
 * \param rng   sample random numbers from this random number generator
 */
char
RNAUtils::sampleNucFromNucDist(const vector<double> &dist, const Runif &rng) {
  double r = rng.runif(0.0, 1.0);
  if (r < dist[RNAUtils::base2int('A')]) return 'A';
  else if (r < dist[RNAUtils::base2int('A')] +\
               dist[RNAUtils::base2int('C')]) return 'C';
  else if (r < dist[RNAUtils::base2int('A')] +\
               dist[RNAUtils::base2int('C')] +\
               dist[RNAUtils::base2int('G')]) return 'G';
  else return 'T';
}

/**
 * \brief sample a single nucleotide from the default background distribution
 * \param rng   sample random numbers from this random number generator
 */
char
RNAUtils::sampleNuc(const Runif &rng) {
  double r = rng.runif(0.0, 1.0);
  if (r < RNAUtils::DEFAULT_BACKGROUND[RNAUtils::base2int('A')])
    return 'A';
  else if (r < RNAUtils::DEFAULT_BACKGROUND[RNAUtils::base2int('A')] +\
               RNAUtils::DEFAULT_BACKGROUND[RNAUtils::base2int('C')])
    return 'C';
  else if (r < RNAUtils::DEFAULT_BACKGROUND[RNAUtils::base2int('A')] +\
               RNAUtils::DEFAULT_BACKGROUND[RNAUtils::base2int('C')] +\
               RNAUtils::DEFAULT_BACKGROUND[RNAUtils::base2int('G')])
    return 'G';
  else return 'T';
}

/**
 * \summary generate a sequence from a nucleotide distribution.
 * \param dist      the distribution to use; dist[base2int('A')] gives the
 *                  probability that an 'A' will be sampled, and so on.
 * \param length    the length of the sequence to build.
 * \param rng       sample random numbers from this random number generator
 * \throw SMITHLABException: if the dist. vector has the wrong dimensions.
 */
string
RNAUtils::sampleSeqFromNucDist(const vector<double> &dist, const size_t length,
                               const Runif &rng) {
  if (dist.size() != RNAUtils::RNA_ALPHABET_SIZE) {
    stringstream ss;
    ss << "Failed to generate sequence from nucleotide distribution, "
       << "distribution vector was malformed: found "
       << dist.size() << " entries; expected " << RNAUtils::RNA_ALPHABET_SIZE;
    throw SMITHLABException(ss.str());
  }

  string res = "";
  for (size_t i = 0; i < length; ++i) res += sampleNucFromNucDist(dist, rng);
  return res;
}

/**
 * \brief           sample a sequence from the default background nucleotide
 *                  distribution.
 * \param length    the length of the sequence to build.
 * \param rng       random number generator to use
 */
string
RNAUtils::sampleSeq(const size_t length, const Runif &rng) {
  string res = "";
  for (size_t i = 0; i < length; ++i) res += RNAUtils::sampleNuc(rng);
  return res;
}

/****
 * \summary     generate a sequence from a position weight matrix
 * \param pwm   the PWM to make the sequence from
 * \param rng   TODO
 * \throw SMITHLABException if the PWM is malformed.
 */
string
RNAUtils::sampleSequenceFromPWM(const vector<vector<double> > pwm,
                                const Runif &rng) {
  string res = "";
  for (size_t j = 0; j < pwm.size(); ++j) {
    if (pwm[j].size() != RNAUtils::RNA_ALPHABET_SIZE) {
      stringstream ss;
      ss << "Failed to generate sequence from PWM, PWM was malformed: found "
         << pwm[j].size() << " entries for position " << j
         << "; expected " << RNAUtils::RNA_ALPHABET_SIZE;
      throw SMITHLABException(ss.str());
    }
    res += sampleNucFromNucDist(pwm[j], rng);
  }
  return res;
}

/**
 * \brief \todo
 * \param seq       \todo
 * \param cnstrnt   \todo
 */
double
RNAUtils::get_minimum_free_energy(const string seq, const string cnstrnt) {

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
void
RNAUtils::get_base_pair_probability_vector(bool VERBOSE,
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
         << "              " << endl;
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
    assert(seqc.size() == conc.size());

    vector<double> bppv_local;
    MC mc;
    mc.init_pf_fold(seqc.size(), energy_pf);
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

