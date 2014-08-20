/*    
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Emad Bahrami-Samani and Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iterator>
#include <numeric>
#include <limits>

#include "smithlab_utils.hpp"
#include "Model.hpp"

using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::accumulate;

// TODO -- PJU there is also an alphabet size specified in RNA_Utils; check
// whether this one is needed and if not, remove.
using smithlab::alphabet_size;

// initialization of non-integral constants
const double Model::pseudocount = 0.1;
const double Model::tolerance = 1e-10;
const double Model::zoops_threshold = 0;
const double Model::DEFAULT_GEO_P = 0.4;

/******************************************************************************
 *        STATIC HELPER FUNCTIONS USED FOR CHECKING DATA CONSISTENCY
 *****************************************************************************/

/***
 * \summary this function checks that the sequences, structure and diagnostic
 *          events vectors are consistent, and throws an exception if they are
 *          not. Consistent means their dimensions are all matching.
 * \param seqs          seqs[i][j] is the nuc. at the jth pos. in the ith seq
 * \param struc         struc[i][j] is the prob. that nuc. j in seq. i is paired
 * \param diagEvents    diagEvents[i][j] is the jth diagnostic event in seq i
 * \param msg           message to use in exception if not consistent
 * \throw SMTIHLABException if the number of structure vectors or the number of
 *                          diagnostic events vectors does not match the number
 *                          of sequences, or if the length of any of those does
 *                          not match the length of the corresponding sequence.
 */
static void
checkAndThrow_consistent(const vector<string> &seqs,
                         const vector<vector<double> > &struc,
                         const vector<vector<double> > &diagEvents,
                         const string &msg) {
  if (diagEvents.size() != seqs.size()) {
      stringstream ss;
      ss << msg << " Expected diagnostic events vector to have same dimension "
         << "as sequences vector (" << seqs.size() << " elements), but it "
         << "doesn't (has " << diagEvents.size() << ")";
      throw SMITHLABException(ss.str());
  }
  if (struc.size() != seqs.size()) {
    stringstream ss;
    ss << msg << " Expected structure vector to have same dimension "
       << "as sequences vector (" << seqs.size() << " elements), but it "
       << "doesn't (has " << struc.size() << ")";
    throw SMITHLABException(ss.str());
  }
  for (size_t i = 0; i < seqs.size(); ++i) {
    if (seqs[i].size() != struc[i].size()) {
      stringstream ss;
      ss << msg << " Expected structure vector for sequence " << i << " to "
         << "have " << seqs[i].size() << " elements, but it doesn't (has "
         << struc[i].size() << ")";
      throw SMITHLABException(ss.str());
    }
  }
}

/******************************************************************************
 *                STATIC HELPER FUNCTIONS USED FOR PRINTING
 ******************************************************************************/

/***
 * \summary return a string representation of a vector of double
 * \param v the vector to turn into a string
 * \TODO -- PJU: should be a template really.
 */
static string
vecToString(const vector<double> &v) {
  stringstream ss;
    for (size_t i = 0; i < v.size(); ++i) {
      ss << v[i];
      if (i != v.size() - 1) ss << ", ";
    }
    return ss.str();
}

/***
 * \summary return a string representation of a matrix of double
 * \param matrix the matrix to turn into a string
 * \TODO -- PJU: should be a template really.
 */
static string
matrixToString(const vector<vector<double> > &matrix) {
  stringstream ss;
  for (size_t i = 0; i < matrix.size(); ++i) {
    ss << vecToString(matrix[i]);
    if (i != matrix.size() - 1) ss << endl;
  }
  return ss.str();
}

/******************************************************************************
 *                 STATIC HELPER FUNCTIONS FOR BASE COUNTING
 ******************************************************************************/

/***
 * \summary produce a probabilistic count of the number of each base that
 *          occurs in the foreground (i.e. in the motif) and background (i.e.
 *          outside of the motif) in a set of sequences, given probabilities on
 *          where the motif occurrences start.
 * \param sequences     sequences[i][j] is the nuc. position j of seq i. The
 *                      sequences may have N's; these won't be counted.
 * \param site_indic    site_indic[i][j] is the prob. that a motif occurs at
 *                      position j of sequence i
 * \param motif_width   the length/width of the motif
 * \param nb_fg         this vector will be populated such that nb_fg[k][b] is
 *                      the probabilistic count of the number of times base b
 *                      is at position k of the motif (0 <= k <= motif_width).
 *                      Any existing entries in the vector will be cleared.
 * \param nb_bg         this vector will be populated such that nb_bg[b]
 *                      is the probabilistic count of the number of times base
 *                      b occurs in the background of the sequences.
 * \pseudoCount         add this pseudo-count to the number of observations for
 *                      each base
 */
static void
calculate_number_of_bases_fg_bg(const vector<string> &seqs,
                                const vector<vector<double> > &site_indic,
                                const size_t motif_width,
                                vector<vector<double> > &nb_fg,
                                vector<double> &nb_bg,
                                const double pseudoCount=1) {
  nb_fg.clear();
  nb_fg.resize(motif_width,
               vector<double>(alphabet_size, pseudoCount));
  for (size_t i = 0; i < site_indic.size(); ++i) {
    for (size_t j = 0; j < site_indic[i].size(); ++j) {
      for (size_t k = 0; k < motif_width; ++k) {
        // can't count anything for N's
        if ((seqs[i][j + k] == 'N') || (seqs[i][j + k] == 'n')) continue;

        const size_t base = base2int(seqs[i][j + k]);
        assert(base < alphabet_size);
        nb_fg[k][base] += site_indic[i][j];
      }
    }
  }

  nb_bg.clear();
  nb_bg.resize(alphabet_size, pseudoCount);
  for (size_t i = 0; i < seqs.size(); ++i) {
    for (size_t j = 0; j < seqs[i].length(); ++j) {
      // can't count anything for N's
      if ((seqs[i][j] == 'N') || (seqs[i][j] == 'n')) continue;

      // motif_width, rather than 1, since each position is considered
      // that many times in the loop below.
      const size_t base = base2int(seqs[i][j]);
      assert(base < alphabet_size);
      nb_bg[base] += motif_width;
    }
  }

  for (size_t i = 0; i < site_indic.size(); ++i) {
    for (size_t j = 0; j < site_indic[i].size(); ++j) {
      for (size_t k = 0; k < motif_width; ++k) {
        // can't count anything for N's
        if ((seqs[i][j + k] == 'N') || (seqs[i][j + k] == 'n')) continue;

        const size_t base = base2int(seqs[i][j + k]);
        assert(base < alphabet_size);
        nb_bg[base] -= site_indic[i][j];
      }
    }
  }
}

/***
* \summary produce a probabilistic count of the number of each base that
 *         occurs in the foreground (i.e. in the motif) and background (i.e.
 *         outside of the motif) with a single-stranded or double-stranded
 *         structure in a set of sequences, given probabilities on
 *         where the motif occurrences start.
 * \param sequences     sequences[i][j] is the nuc. position j of seq i. The
 *                      sequences may have N's; these won't be counted.
 * \param secStr        TODO
 * \param site_indic    site_indic[i][j] is the prob. that a motif occurs at
 *                      position j of sequence i
 * \param motif_width   the length/width of the motif
 * \param nb_fg_ss      TODO
 * \param nb_bg_ss      TODO
 * \param nb_fg_ds      TODO
 * \param nb_bg_ds      TODO
 * \pseudoCount         add this pseudo-count to the number of observations for
 *                      each base
 */
static void
calculate_number_of_bases_fg_bg_str(const vector<string> &seqs,
                                    const vector<vector<double> > &secStr,
                                    const vector<vector<double> > &site_indic,
                                    const size_t motif_width,
                                    vector<vector<double> > &nb_fg_ss,
                                    vector<double> &nb_bg_ss,
                                    vector<vector<double> > &nb_fg_ds,
                                    vector<double> &nb_bg_ds,
                                    const double pseudoCount=1) {

  nb_fg_ss.clear();
  nb_fg_ss.resize(motif_width, vector<double>(alphabet_size, pseudoCount));
  nb_fg_ds.clear();
  nb_fg_ds.resize(motif_width, vector<double>(alphabet_size, pseudoCount));

  // this loop calculates the number of bases occurring at a
  // particular location of the motif with particular secondary
  // structure
  for (size_t i = 0; i < site_indic.size(); ++i) {
    for (size_t j = 0; j < site_indic[i].size(); ++j) {
      for (size_t k = 0; k < motif_width; ++k) {
        // can't count anything for N's
        if ((seqs[i][j + k] == 'N') || (seqs[i][j + k] == 'n')) continue;

        const size_t base_idx = base2int(seqs[i][j + k]);
        assert(base_idx < alphabet_size);

        const double curr_sec_str = secStr[i][j + k];
        const double curr_site_indic = site_indic[i][j];
        nb_fg_ss[k][base_idx] += (1.0 - curr_sec_str) * curr_site_indic;
        nb_fg_ds[k][base_idx] += curr_sec_str * curr_site_indic;
      }
    }
  }

  nb_bg_ss.clear();
  nb_bg_ss.resize(alphabet_size, pseudoCount);
  nb_bg_ds.clear();
  nb_bg_ds.resize(alphabet_size, pseudoCount);

  // the next two loop sets calculate the number of bases occurring in
  // the background with a particular secondary structure
  for (size_t i = 0; i < seqs.size(); ++i) {
    for (size_t j = 0; j < seqs[i].length(); ++j) {
      // can't count anything for N's
      if ((seqs[i][j] == 'N') || (seqs[i][j] == 'n')) continue;

      const size_t base_idx = base2int(seqs[i][j]);
      assert(base_idx < alphabet_size);
      nb_bg_ss[base_idx] += (1.0 - secStr[i][j]) * motif_width;
      nb_bg_ds[base_idx] += (secStr[i][j] * motif_width);
    }
  }
  for (size_t i = 0; i < site_indic.size(); ++i) {
    for (size_t j = 0; j < site_indic[i].size(); ++j) {
      for (size_t k = 0; k < motif_width; ++k) {
        // can't count anything for N's
        if ((seqs[i][j + k] == 'N') || (seqs[i][j + k] == 'n')) continue;

        const size_t base_idx = base2int(seqs[i][j + k]);
        assert(base_idx < alphabet_size);
        nb_bg_ss[base_idx] -= ((1.0 - secStr[i][j + k]) * site_indic[i][j]);
        nb_bg_ds[base_idx] -= (secStr[i][j + k] * site_indic[i][j]);
      }
    }
  }
}


/******************************************************************************
 *          STATIC HELPER FUNCTIONS USED DURING MODEL MAXIMIZATION
 ******************************************************************************/

/***
 * \summary maximize the sequence component of the model given a set of
 *          sequences and indicator variables specifying probabilities of where
 *          the motif occurs in them.
 * \param seqs      TODO
 * \param siteInd   TODO
 * \param seqIndic  TODO
 * \param matrix    TODO
 * \param freq      TODO
 * \param gamma     TODO
 */
static void
maximization_seq(const vector<string> &seqs,
                 const vector<vector<double> > &siteInd,
                 vector<double> &seqIndic,
                 vector<vector<double> > &matrix,
                 vector<double> &freq,
                 double &gamma) {
  static const double TINY = 1e-100;

  vector<vector<double> > nb_fg;
  vector<double> nb_bg;
  calculate_number_of_bases_fg_bg(seqs, siteInd, matrix.size(), nb_fg, nb_bg);

  for (size_t i = 0; i < matrix.size(); ++i) {
    const double total = accumulate(nb_fg[i].begin(), nb_fg[i].end(), 0.0);
    assert(std::isfinite(total));
    transform(nb_fg[i].begin(), nb_fg[i].end(), matrix[i].begin(),
              std::bind2nd(std::divides<double>(), total));

    // don't let anything get to zero; that'll be bad for log later.
    // also, check that everything is finite.
    for (size_t j = 0; j < matrix[i].size(); ++j) {
      if (matrix[i][j] < TINY) matrix[i][j] = TINY;
      assert(std::isfinite(matrix[i][j]));
    }
  }

  const double total = accumulate(nb_bg.begin(), nb_bg.end(), 0.0);
  transform(nb_bg.begin(), nb_bg.end(), freq.begin(),
            std::bind2nd(std::divides<double>(), total));

  // set gamma, which is, in essence, the probability that all sequences
  // in the data have an occurrence of the motif, assuming all prior probs of
  // seeing a motif in any sequence are equal. Be careful here to not let it
  // exceed 1.0, or get to exactly 0 (we take log of it later)
  gamma = accumulate(seqIndic.begin(), seqIndic.end(), 0.0) / seqs.size();
  gamma = std::max(std::min(gamma, 1.0), TINY);
}

/***
 * \summary maximize parameters for the structure component of the model based
 *          on the data and the indicators.
 * \param secondary_structure   secondary_structure[i][j] gives the probability
 *                              that the jth nuc. in the ith sequence is paired
 * \param site_indic            site_indic[i][j] gives the probability that the
 *                              jth position within the ith sequence contains
 *                              and occurrence of the motif
 * \param seq_indic             TODO
 * \param matrix                TODO
 * \param motif_sec_str         motif_sec_str[i] gives the probability that the
 *                              ith position within the motif is double-stranded.
 * \param f_sec_str             the background base-pair probability
 */
static void
maximization_str(const vector<vector<double> > &secondary_structure,
                 const vector<vector<double> > &site_indic,
                 const vector<double> &seq_indic,
                 vector<vector<double> > &matrix,
                 vector<double> &motif_sec_str,
                 double &f_sec_str) {

  motif_sec_str.clear();
  motif_sec_str.resize(matrix.size(), 0.0);
  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < site_indic.size(); ++j) {
      for (size_t site = 0; site < site_indic[j].size(); ++site) {
        motif_sec_str[i] += seq_indic[j] * site_indic[j][site]
                            * secondary_structure[j][site + i];
      }
    }

    motif_sec_str[i] = motif_sec_str[i]
        / accumulate(seq_indic.begin(), seq_indic.end(), 0.0);
  }

  // if we want to learn this parameter we have to calculate this here, but for
  // the moment we just set it to be always 0.5, based on empirical observation
  f_sec_str = 0.5;
}

static double 
newtonRaphson (const vector<vector<double> > &diagEvents,
               const vector<vector<double> > &siteInd,
               const double &geoP,
               const int &geoDelta)
{
  double numerator = 0.0;
  double denominator = 0.0;
  for (size_t i = 0; i < siteInd.size(); i++)
    for (size_t j = 0; j < siteInd[i].size(); j++) {
      vector<double> sum_1;
      vector<double> sum_2;
      vector<double> sum_3;
      for (size_t l = 0; l < diagEvents[i].size(); l++) {
        double term_1 = 0.0;
        double term_2 = 0.0;
        double term_3 = 0.0;
        double term_4 = 0.0;
        double term_5 = 0.0;

        if (abs(l - (j + geoDelta)) > 0)
          term_1 = log(abs(l - (j + geoDelta)));

        if ((abs(l - (j + geoDelta)) - 1) > 0)
          term_2 = log(abs(l - (j + geoDelta)) - 1.0);

        if ((abs(l - (j + geoDelta)) - 2) > 0)
          term_3 = ((abs(l - (j + geoDelta)) - 2) * log(1.0 - geoP));

        if (abs(l - (j + geoDelta)) > 0)
          term_4 = (abs(l - (j + geoDelta)) * log(1.0 - geoP));

        if ((abs(l - (j + geoDelta)) - 1) > 0)
          term_5 = ((abs(l - (j + geoDelta)) - 1) * log(1.0 - geoP));

        sum_1.push_back(log(diagEvents[i][l]) + term_1 + term_2 + term_3);
        sum_2.push_back(log(diagEvents[i][l]) + term_4);
        sum_3.push_back(log(diagEvents[i][l]) + term_1 + term_5);
      }
      const double log_A = smithlab::log_sum_log_vec(sum_1, sum_1.size());
      const double log_B = smithlab::log_sum_log_vec(sum_2, sum_2.size());
      const double log_C = smithlab::log_sum_log_vec(sum_3, sum_3.size());
      const double num_ij = (1.0/geoP) - (exp(log_C - log_B));
      const double denom_ij = (-1.0/(geoP*geoP) - (exp(log_A - log_B) - exp(2.0*(log_C - log_B))));
      numerator += (siteInd[i][j] * num_ij);
      denominator += (siteInd[i][j] * denom_ij);
    }
  return geoP - (numerator / denominator);
}

/***
 * \summary given a set of diagnostic event locations within sequences, and
 *          both indicators on the probability of motif occurrence at each
 *          position within each sequence and the probability that each
 *          sequence has a motif, calculate the MLE estimate for the geometric
 *          distribution of distances between diagnostic events and motif
 *          occurrences.
 * \param diagEvents diagEvents[i][k] is the (relative) position of the kth
 *                   diagnostic event in sequence i
 * \param siteInd    siteInd[i][j] is the probability that the motif occurrence
 *                   in sequence i is at position j, given that an sequence i
 *                   contains an occurrence.
 * \param seqInd     seqInd[i] is the prob. that sequence i contains an
 *                   occurrence of the motif
 * \param matrix     matrix[j][b] is the prob. that base b occurs at position j
 *                   of the motif
 * \param geoP       TODO
 * \param geoDelta   TODO
 * \bug \todo PJU: geoDelta isn't updated here -- should it be? if not, it
 *            should be set to const, and an explanation about why it isn't
 *            updated added to comments.
 */
static void
maximization_de(const vector<vector<double> > &diagEvents,
                const vector<vector<double> > &siteInd,
                const vector<double> &seqInd,
                vector<vector<double> > &matrix,
                double &geoP,
                int &geoDelta) {

  geoP = Model::DEFAULT_GEO_P;
  const double tol = 0.0001;    // 0.0001 is the error level we wish
  double old;

  do
  {
    old = geoP;
    geoP = max(min(newtonRaphson(diagEvents, siteInd, geoP, geoDelta), 0.999),
               std::numeric_limits<double>::min());
  }
  while (abs(old - geoP) > tol);     // while loop because the

//  if (!std::isfinite(geoP)) {
//    stringstream ss;
//    ss << "failed maximization of geometric parameter; numerator was "
//       << matrixToString(siteInd) << " seqInd were " << vecToString(seqInd);
//    throw SMITHLABException(ss.str());
//  }
}


/******************************************************************************
 *   STATIC HELPER FUNCTIONS USED IN CALCULATION OF INDICATORS / EXPECTATION
 ******************************************************************************/

/***
 * \summary TODO
 */
static void
get_numerator_seq_de_for_site(const string &seq,
                              const vector<double> &diagEvents,
                              const vector<vector<double> > &matrix,
                              const vector<double> &freqs,
                              const double geo_p,
                              const int geo_delta,
                              const double gamma,
                              const size_t site,
                              double &num) {
  // the log-likelihood of a nucleotide in the foreground (motif) if it is
  // an N. We make this pretty bad, but be careful to make sure we could sum
  // this across all positions of the matrix without causing an overflow.
  const double N_LOG_PROB = -10000;

  vector<double> f_powers(alphabet_size, 0.0);
  for (size_t i = 0; i < seq.length(); ++i) {
    const size_t base = base2int(seq[i]);
    if (i >= site && i < site + matrix.size()) {
      // Penalize N's heavily if they fall in the motif
      if ((seq[i] == 'N') || (seq[i] == 'n')) {
        num += N_LOG_PROB;
      }
      else {
        assert(base < alphabet_size);
        num += log(matrix[i - site][base]);
        if (!std::isfinite(num)) {
          stringstream ss;
          ss << "failed expectation calculation; numerator non-finite. Matrix "
             << "entry was: " << matrix[i - site][base];
          throw SMITHLABException(ss.str());
        }
      }
    } else {
      // Just ignore N's outside of the motif
      if ((seq[i] != 'N') && (seq[i] != 'n')) {
        assert(base < alphabet_size);
        f_powers[base]++;
      }
      if (!std::isfinite(f_powers[base]))
        throw SMITHLABException("failed expectation calc; f_powers non-finite");
    }
  }

  for (size_t b = 0; b < alphabet_size; b++) {
    num += f_powers[b] * log(freqs[b]);
    assert(std::isfinite(num));
  }

  if (diagEvents.size() > 0) { 
    vector<double> powers;
    for (size_t j = 0; j < seq.length(); j++)
      powers.push_back(log(diagEvents[j]) + log(geo_p) + (abs(j - (site + geo_delta)) * log(1.0-geo_p)));
    num += smithlab::log_sum_log_vec(powers, powers.size());
  }
  num += log(gamma);
}

/***
 * \summary TODO
 */
static void
get_numerator_seq_str_de_for_site(const string &seq,
                               const vector<double> &secondary_structure,
                               const vector<double> &diagnostic_events,
                               const vector<vector<double> > &matrix,
                               const vector<double> &motif_sec_str,
                               const vector<double> &freqs,
                               const double f_sec_str,
                               const double geo_p,
                               const int geo_delta,
                               const double gamma,
                               const size_t site,
                               double &num) {
  const double N_LOG_PROB = -10000;

  // calculating the contribution of the foreground and the powers
  // that will be needed for the background calculations (below).
  vector<double> f_powers_ss(alphabet_size, 0.0);
  vector<double> f_powers_ds(alphabet_size, 0.0);
  for (size_t i = 0; i < seq.length(); ++i) {
    const size_t base = base2int(seq[i]);
    if (i >= site && i < site + matrix.size()) {
      // Penalize N's heavily if they fall in the motif
      if ((seq[i] == 'N') || (seq[i] == 'n')) {
        num += N_LOG_PROB;
      }
      else {
        assert(base < alphabet_size);
        num += (secondary_structure[i]
            * log(matrix[i - site][base] * motif_sec_str[i - site]));
        num += ((1.0 - secondary_structure[i])
            * log(matrix[i - site][base] * (1.0 - motif_sec_str[i - site])));
        if (!std::isfinite(num)) {
          stringstream ss;
          ss << "failed expectation calculation; numerator non-finite. BPP was: "
             << secondary_structure[i] << " matrix entry was: "
             << matrix[i - site][base] << " motif secondary structure was: "
             << motif_sec_str[i - site];
          throw SMITHLABException(ss.str());
        }
      }
    } else {
      // Just ignore N's outside of the motif
      if ((seq[i] != 'N') && (seq[i] != 'n')) {
        assert(base < alphabet_size);
        f_powers_ss[base] += (1.0 - secondary_structure[i]);
        f_powers_ds[base] += (secondary_structure[i]);
      }
    }
    if (!std::isfinite(f_powers_ss[base]))
      throw SMITHLABException("failed expectation calc; f_powers_ss non-finite");
    if (!std::isfinite(f_powers_ds[base]))
      throw SMITHLABException("failed expectation calc; f_powers_ds non-finite");
  }

  // calculating the contribution of the background (outside the motif
  // occurrences)
  for (size_t b = 0; b < alphabet_size; b++) {
    num += f_powers_ss[b] * log(freqs[b] * (1.0 - f_sec_str));
    num += f_powers_ds[b] * log(freqs[b] * (f_sec_str));
    assert(std::isfinite(num));
  }

  if (diagnostic_events.size() > 0) {
    vector<double> powers;
    for (size_t j = 0; j < seq.length(); j++) 
      powers.push_back(log(diagnostic_events[j]) + log(geo_p) + (abs(j - (site + geo_delta)) * log(1.0-geo_p)));
    num += smithlab::log_sum_log_vec(powers, powers.size());
  }

  num += log(gamma);
}

/***
 * \summary TODO
 * \param gamma the fraction of the sequences that contain the motif
 */
static void
expectation_seq_str_de_for_single_seq(const string &seq,
                                      const vector<double> &secondary_structure,
                                      const vector<double> &diagnostic_events,
                                      const vector<vector<double> > &matrix,
                                      const vector<double> &motif_sec_str,
                                      const vector<double> &freqs,
                                      const double f_sec_str,
                                      const double geo_p,
                                      const int geo_delta,
                                      const double gamma,
                                      vector<double> &site_indic,
                                      double &seq_indic) {
  // used to stop values reaching exactly zero and then taking their log..
  const double TINY = 1e-100;

  // sanity check on gamma
  if ((gamma > 1.0 ) || (gamma < 0.0)) {
    stringstream ss;
    ss << "failed expectation step: gamma (faction of sequences containing "
       << "the motif) was outside the expected bounds: " << gamma;
    throw SMITHLABException(ss.str());
  }

  // get log likelihood for each site
  vector<double> numerator(site_indic.size(), 0.0);
  for (size_t i = 0; i < site_indic.size(); ++i) {
    get_numerator_seq_str_de_for_site(seq, secondary_structure,
                                      diagnostic_events, matrix, motif_sec_str,
                                      freqs, f_sec_str, geo_p, geo_delta, gamma,
                                      i, numerator[i]);
    assert(std::isfinite(numerator[i]));
  }

  double no_motif = 0.0;
  for (size_t i = 0; i < seq.length(); i++) {
    // skip N's
    if ((seq[i] == 'N') || (seq[i] == 'n')) continue;
    const size_t base = base2int(seq[i]);
    assert(base < alphabet_size);

    no_motif += secondary_structure[i]
        * log(freqs[base] * f_sec_str);
    no_motif += ((1.0 - secondary_structure[i])
        * log(freqs[base] * (1.0 - f_sec_str)));
  }

/*  if (diagnostic_events.size() > 0) {
    vector<double> powers;
    for (size_t site = 0; site < site_indic.size(); ++site)
      for (size_t j = 0; j < seq.length(); j++)
        powers.push_back(log((1.0/site_indic.size())) + log(diagnostic_events[j]) + log(geo_p) + (abs(j - (site + geo_delta)) * log(1.0-geo_p)));
    no_motif += smithlab::log_sum_log_vec(powers, powers.size());;
  }*/

  const double fracSeqsWithoutMotif = std::min((1.0 - gamma) + TINY, 1.0);
  numerator.push_back(no_motif + log(fracSeqsWithoutMotif));

  const double denominator = smithlab::log_sum_log_vec(numerator,
                                                       numerator.size());
  assert(std::isfinite(denominator));
  for (size_t i = 0; i < site_indic.size(); ++i)
    site_indic[i] = exp(numerator[i] - denominator);

  seq_indic = accumulate(site_indic.begin(), site_indic.end(), 0.0);
  assert(std::isfinite(seq_indic));
}

/***
 * \summary TODO
 * \param gamma the fraction of the sequences that contain the motif
 */
static void
expectation_seq_de_for_single_seq(const string &seq,
                                  const vector<double> &diagnostic_events,
                                  const vector<vector<double> > &matrix,
                                  const vector<double> &freqs,
                                  const double geo_p,
                                  const int geo_delta,
                                  const double gamma,
                                  vector<double> &site_indic,
                                  double &seq_indic) {
  // used to stop values reaching exactly zero and then taking their log..
  const double TINY = 1e-100;

  // sanity check on gamma
  if ((gamma > 1.0 ) || (gamma < 0.0)) {
    stringstream ss;
    ss << "failed expectation step: gamma (faction of sequences containing "
       << "the motif) was outside the expected bounds: " << gamma;
    throw SMITHLABException(ss.str());
  }

  // get log likelihood for each site
  vector<double> numerator(site_indic.size(), 0.0);
  for (size_t i = 0; i < site_indic.size(); ++i) {
    get_numerator_seq_de_for_site(seq, diagnostic_events, matrix, freqs, geo_p,
                                  geo_delta, gamma, i, numerator[i]);
    assert(std::isfinite(numerator[i]));
  }

  double no_motif = 0.0;
  for (size_t i = 0; i < seq.length(); i++) {
    // ignore N's
    if ((seq[i] == 'N') || (seq[i] == 'n')) continue;
    const size_t base = base2int(seq[i]);
    assert(base < alphabet_size);
    no_motif += log(freqs[base]);
  }

/*   if (diagnostic_events.size() > 0) {
    vector<double> powers;
    for (size_t site = 0; site < site_indic.size(); ++site)
      for (size_t j = 0; j < seq.length(); j++)
        powers.push_back(log((1.0/site_indic.size())) + log(diagnostic_events[j]) + log(geo_p) + (abs(j - (site + geo_delta)) * log(1.0-geo_p)));
    no_motif += smithlab::log_sum_log_vec(powers, powers.size());
  }*/

  const double fracSeqsWithoutMotif = std::min((1.0 - gamma) + TINY, 1.0);
  numerator.push_back(no_motif + log(fracSeqsWithoutMotif));

  const double denominator = smithlab::log_sum_log_vec(numerator,
                                                       numerator.size());
  for (size_t i = 0; i < site_indic.size(); ++i) {
    site_indic[i] = exp(numerator[i] - denominator);
    assert(std::isfinite(site_indic[i]));
  }

  seq_indic = accumulate(site_indic.begin(), site_indic.end(), 0.0);
}

/***
 * \summary calculate the occurrence indicators for a sequence given a model
 *          of the motif. This function considers only the sequence, not the
 *          structure nor any diagnostic events.
 * \param seq           TODO
 * \param matrix        TODO
 * \param freqs         TODO
 * \param gamma         TODO
 * \param site_indic    TODO
 * \param seq_indic     TODO
 */
static void
expectation_for_single_seq(const string &seq,
                           const vector<vector<double> > &matrix,
                           const vector<double> &freqs,
                           const double gamma,
                           vector<double> &site_indic,
                           double &seq_indic) {
  // sequence expectation is equivalent to sequence and DE expectation, but
  // with no DEs, so...
  vector<double> diagEvents;
  // geoP and geoDelta could be anything really, they'll be passed through,
  // but won't be used since diagEvents.size() == 0
  const double geoP = 1;
  const int geoDelta = 0;
  expectation_seq_de_for_single_seq(seq, diagEvents, matrix, freqs, geoP,
                                    geoDelta, gamma, site_indic, seq_indic);
}

/***
 * \summary TODO
 */
static void
expectation_seq(const vector<string> &sequences,
                const vector<vector<double> > &matrix,
                const vector<double> &freqs,
                const double gamma,
                vector<vector<double> > &site_indic,
                vector<double> &seq_indic) {

  for (size_t i = 0; i < sequences.size(); i++)
    expectation_for_single_seq(sequences[i], matrix, freqs, gamma,
                               site_indic[i], seq_indic[i]);
}

/***
 * \summary TODO
 */
static void
expectation_seq_de(const vector<string> &sequences,
                   const vector<vector<double> > &diagnostic_events,
                   const vector<vector<double> > &matrix,
                   const vector<double> &freqs,
                   const double geo_p,
                   const int geo_delta,
                   const double gamma,
                   vector<vector<double> > &site_indic,
                   vector<double> &seq_indic) {
  if (Model::DEBUG_LEVEL >= 2)
    cerr << "performing expectation step with matrix " << endl
         << matrixToString(matrix) << endl;
  for (size_t i = 0; i < sequences.size(); i++)
    expectation_seq_de_for_single_seq(sequences[i], diagnostic_events[i],
                                      matrix, freqs, geo_p, geo_delta, gamma,
                                      site_indic[i], seq_indic[i]);
  if (Model::DEBUG_LEVEL >= 2)
    cerr << "finished expectation step" << endl;
}

/***
 * \summary TODO
 */
static void
expectation_seq_str_de(const vector<string> &sequences,
                       const vector<vector<double> > &secondary_structure,
                       const vector<vector<double> > &diagnostic_events,
                       const vector<vector<double> > &matrix,
                       const vector<double> &motif_sec_str,
                       const vector<double> &freqs,
                       const double f_sec_str,
                       const double geo_p,
                       const int geo_delta,
                       const double gamma,
                       vector<vector<double> > &site_indic,
                       vector<double> &seq_indic) {
  for (size_t i = 0; i < sequences.size(); i++)
    expectation_seq_str_de_for_single_seq(sequences[i], secondary_structure[i],
                                          diagnostic_events[i], matrix,
                                          motif_sec_str, freqs, f_sec_str,
                                          geo_p, geo_delta, gamma,
                                          site_indic[i], seq_indic[i]);
}

/******************************************************************************
 *                             SIMPLE GETTERS
 *****************************************************************************/
string
Model::toString_pwm() const {
  return matrixToString(this->matrix);
}

/******************************************************************************
 *                GETTERS -- LIKELIHOOD CALCULATION FUNCTIONS
 ******************************************************************************/

/***
 * \summary TODO
 */
double
Model::calculate_oops_log_l(const vector<string> &sequences,
                            const vector<vector<double> > &site_indic) const {

  vector<vector<double> > nb_fg;
  vector<double> nb_bg;
  calculate_number_of_bases_fg_bg(sequences, site_indic, matrix.size(), nb_fg,
                                  nb_bg);

  double ret = 0.0;
  for (size_t i = 0; i < alphabet_size; ++i) {
    ret += nb_bg[i] * log(f[i]);
    for (size_t j = 0; j < matrix.size(); ++j)
      ret += nb_fg[j][i] * log(matrix[j][i]);
  }

  for (size_t i = 0; i < site_indic.size(); ++i)
    ret -= log(site_indic[i].size());

  return ret;
}

/***
 * \summary Calculate the zoops log likelihood for this model given a set of
 *          sequences and indicators for motif occurrences in those sequences.
 *          This function ignores secondary structure and diagnostic events.
 * \param seqs          TODO
 * \param site_indic    site_indic[i][j] is the probability that the
 *                      jth position in the ith sequence is the start
 *                      of a motif occurrence.
 * \param seq_indic     TODO
 */
double
Model::calculate_zoops_log_l(const vector<string> &seqs,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {
  // calc. without diag. events is equivalent to having no diag. events, so...
  vector<vector<double> > diagEvents (seqs.size());
  return calculate_zoops_log_l(seqs, diagEvents, site_indic, seq_indic);
}

/***
 * \summary Calculate the zoops log-likelihood for this model given a set of
 *          sequences, diganostic events that occur in those sequences, and
 *          indicators for occurrences of the motif in the sequences.
 *          This function ignores secondary structure.
 * \param sequences             TODO
 * \param diagnostic_events     TODO
 * \param site_indic            site_indic[i][j] is the probability that the
 *                              jth position in the ith sequence is the start
 *                              of a motif occurrence.
 * \param seq_indic             TODO
 */
double
Model::calculate_zoops_log_l(const vector<string> &sequences,
                             const vector<vector<double> > &diagnostic_events,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {
  // used to stop values reaching exactly zero and then taking their log..
  const double TINY = 1e-100;

  vector<vector<double> > nb_fg;
  vector<double> nb_bg;
  calculate_number_of_bases_fg_bg(sequences, site_indic, matrix.size(), nb_fg,
                                  nb_bg);

  double ret = 0.0;
  for (size_t i = 0; i < alphabet_size; ++i) {
    ret += nb_bg[i] * log(f[i]);
    if (!std::isfinite(ret)) {
      stringstream ss;
      ss << "log-likelihood calculation failed; result not finite. background "
         << "count for base was " << nb_bg[i] << "; f for base was " << f[i];
      throw SMITHLABException(ss.str());
    }
    for (size_t j = 0; j < matrix.size(); ++j) {
      ret += nb_fg[j][i] * log(matrix[j][i]);
      if (!std::isfinite(ret)) {
        stringstream ss;
        ss << "log-likelihood calculation failed; result not finite. "
           << "foreground count for base was " << nb_fg[j][i] << "; matrix "
           << "entry for base was " << matrix[j][i];
        throw SMITHLABException(ss.str());
      }
    }
  }
  if (diagnostic_events.size() > 0) {
    for (size_t i = 0; i < sequences.size(); i++) {
      if (diagnostic_events[i].size() > 0) {
        for (size_t k = 0; k < site_indic[i].size(); k++) {
          vector<double> powers;
          for (size_t j = 0; j < sequences[i].length(); j++)
            powers.push_back(log(diagnostic_events[i][j]) + log(p) + (abs(j - (k + delta)) * log(1.0-p)));
          ret += (site_indic[i][k] * smithlab::log_sum_log_vec(powers, powers.size()));
        }
      }
    }
  }

  //--------
  for (size_t i = 0; i < sequences.size(); i++) {
    double has_no_motif = 0.0;
    for (size_t j = 0; j < sequences[i].length(); j++) {
      // ignore N's
      if ((sequences[i][j] == 'N') || (sequences[i][j] == 'n')) continue;
      const size_t base = base2int(sequences[i][j]);
      assert(base < alphabet_size);
      has_no_motif += log(f[base]);
    }

/*    if (diagnostic_events[i].size() > 0) {
      vector<double> powers;
      for (size_t k = 0; k < site_indic[i].size(); k++) {
        for (size_t j = 0; j < sequences[i].length(); j++) {
          powers.push_back(log((1.0/site_indic.size())) + log(diagnostic_events[i][j]) + log(p) + (abs(j - (k + delta)) * log(1.0-p)));
        }
      }
      has_no_motif += smithlab::log_sum_log_vec(powers, powers.size());
    }
*/
    ret += (1.0 - seq_indic[i]) * has_no_motif;
    ret += (1.0 - seq_indic[i]) * log(std::min(1.0 - gamma + TINY, 1.0));
    ret += seq_indic[i] * log(gamma);
  }
  //--------

  return ret;
}

/***
 * \summary calculate the zoops log-likelihood for this model given a set of
 *          sequences, secondary structure and diagnostic events.
 */
double
Model::calculate_zoops_log_l(const vector<string> &sequences,
                             const vector<vector<double> > &secondary_structure,
                             const vector<vector<double> > &diagnostic_events,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {
  // used to stop values reaching exactly zero and then taking their log..
  const double TINY = 1e-100;

  // this function makes some implicit assumptions; let's just check them..
  checkAndThrow_consistent(sequences, secondary_structure,
                           diagnostic_events, "Likelihood calculation failed.");

  vector<vector<double> > nb_fg_ss;
  vector<vector<double> > nb_fg_ds;
  vector<double> nb_bg_ss;
  vector<double> nb_bg_ds;
  calculate_number_of_bases_fg_bg_str(sequences, secondary_structure,
                                      site_indic, matrix.size(), nb_fg_ss,
                                      nb_bg_ss, nb_fg_ds, nb_bg_ds);

  // this loop calculates the log likelihood of the model using the
  // counts obtained from the previous function call. This loop is the
  // same as would work in an OOPS model, ZOOPS-specific stuff is
  // below.
  double ret = 0.0;
  for (size_t i = 0; i < alphabet_size; ++i) {
    ret += (nb_bg_ss[i] * log(f[i] * (1.0 - f_sec_str)));
    ret += (nb_bg_ds[i] * log(f[i] * f_sec_str));
    for (size_t j = 0; j < matrix.size(); ++j) {
      ret += (nb_fg_ss[j][i] * log(matrix[j][i] * (1.0 - motif_sec_str[j])));
      ret += (nb_fg_ds[j][i] * log(matrix[j][i] * motif_sec_str[j]));
    }
  }

  // ZOOPS-specific calculation starts from here --------
  if (diagnostic_events.size() > 0) {
    for (size_t i = 0; i < sequences.size(); i++) {
      if (diagnostic_events[i].size() > 0) {
        for (size_t k = 0; k < site_indic[i].size(); k++) {
          vector<double> powers;
          for (size_t j = 0; j < sequences[i].length(); j++)
            powers.push_back(log(diagnostic_events[i][j]) + log(p) + (abs(j - (k + delta)) * log(1.0-p)));
          ret += (site_indic[i][k] * smithlab::log_sum_log_vec(powers, powers.size()));
        }
      }
    }
  }

  for (size_t i = 0; i < sequences.size(); i++) {
    double has_no_motif = 0.0;
    for (size_t j = 0; j < sequences[i].length(); j++) {
      // ignore N's
      if ((sequences[i][j] == 'N') || (sequences[i][j] == 'n')) continue;
      const size_t base = base2int(sequences[i][j]);
      assert(base < alphabet_size);
      has_no_motif += log(f[base]);
    }

/*    if (diagnostic_events[i].size() > 0) {
      vector<double> powers;
      for (size_t k = 0; k < site_indic[i].size(); k++) {
        for (size_t j = 0; j < sequences[i].length(); j++) {
          powers.push_back(log((1.0/site_indic.size())) + log(diagnostic_events[i][j]) + log(p) + (abs(j - (k + delta)) * log(1.0-p)));
        }
      }
      has_no_motif += smithlab::log_sum_log_vec(powers, powers.size());
    }*/

    ret += (1.0 - seq_indic[i]) * has_no_motif;
    ret += (1.0 - seq_indic[i]) * log(std::min(1.0 - gamma + TINY, 1.0));
    ret += seq_indic[i] * log(gamma);
  }
  // ZOOPS-specific calculation ends here --------

  return ret;
}


/******************************************************************************
 * SETTERS -- PARAMETER ESTIMATION BY EXPECTATION MAXIMIZATION
 ******************************************************************************/

/***
 * \summary use the expectation maximization algorithm to estimate the
 *          parameters for this model given a set of sequences and their
 *          structure
 * \param sequences         TODO
 * \param secStr            TODO
 * \param site_indic        TODO
 * \param seq_indic         TODO
 */
void
Model::expectation_maximization_seq_str(const vector<string> &sequences,
                                        const vector<vector<double> > &secStr,
                                        vector<vector<double> > &site_indic,
                                        vector<double> &seq_indic) {
  // performing EM with sequence and structure but not diagnostic events is
  // equivalent to doing it with diagnsotic event vectors, but with them all
  // empty, so..
  vector<vector<double> > diagEvents(sequences.size());
  expectationMax_SeqStrDE(sequences, secStr, diagEvents,
                          site_indic, seq_indic);
}

/***
 * \summary use the expectation maximization algorithm to estimate the
 *          parameters for this model given a set of sequences
 */
void
Model::expectation_maximization_seq(const vector<string> &seqs,
                                    vector<vector<double> > &site_indic,
                                    vector<double> &seq_indic) {
  bool first = true;
  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq(seqs, matrix, f, gamma, site_indic, seq_indic);
    maximization_seq(seqs, site_indic, seq_indic, matrix, f, gamma);
    const double score = calculate_zoops_log_l(seqs, site_indic, seq_indic);
    if (!first) {
      const double delta = fabs(prev_score - score);
      const double deltaProp = delta / fabs(prev_score);
      if (Model::DEBUG_LEVEL >= 2)
        cerr << "new ll: " << score << " old ll: " << prev_score << " "
             << "absolute change is " << delta << " fraction change "
             << deltaProp << endl;
      if (deltaProp < tolerance) break;
    } else {
      first = false;
    }
    prev_score = score;
  }
}

/***
 * \summary use the expectation maximization algorithm to estimate the
 *          parameters for this model given a set of sequences, and diagnostic
 *          events within them.
 * \param seqs         TODO
 * \param diagEvents   diagnostic_events[i][j] is the location of the jth
 *                     diagnostic event in the ith sequence.
 *                     diagnostic_events[i] can be empty for any value of i
 *                     (i.e. there are no diagnostic events in that sequence)
 * \param site_indic   TODO
 * \param seq_indic    TODO
 * \param holdDelta    if true, delta (the offset of the geometric mode from
 *                     the start of the motif) is not updated. This parameter
 *                     has a default value of false if not provided.
 */
void
Model::expectation_maximization_seq_de(const vector<string> &seqs,
                                       const vector<vector<double> > &diagEvents,
                                       vector<vector<double> > &site_indic,
                                       vector<double> &seq_indic,
                                       const bool holdDelta = false) {
  if (!holdDelta) this->estimateDelta(seqs, diagEvents);
  double prev_score = std::numeric_limits<double>::max();
  bool first = true;
  double score = 0.0;
  for (size_t i = 0; i < max_iterations; ++i) {
    if (Model::DEBUG_LEVEL >= 1) {
      cerr << "EM, sequence and DE, iteration number " << i << endl
           << "expectation step" << endl;
    }
    expectation_seq_de(seqs, diagEvents, matrix, f, p, delta, gamma,
                       site_indic, seq_indic);
    if (Model::DEBUG_LEVEL >= 1) cerr << "sequence maximization step" << endl;
    maximization_seq(seqs, site_indic, seq_indic, matrix, f, gamma);
    if (Model::DEBUG_LEVEL >= 1) cerr << "DE maximization step" << endl;
    maximization_de(diagEvents, site_indic, seq_indic, matrix, p, delta);
    if (Model::DEBUG_LEVEL >= 1) cerr << "calculating log-likelihood" << endl;
    score = calculate_zoops_log_l(seqs, diagEvents, site_indic, seq_indic);
    if (!first) {
      const double delta = fabs(prev_score - score);
      const double deltaProp = delta / fabs(prev_score);
      if (Model::DEBUG_LEVEL >= 1)
        cerr << "new ll: " << score << " old ll: " << prev_score << " "
             << "absolute change is " << delta << " fraction change "
             << deltaProp << endl;
      if (deltaProp < tolerance) break;
    } else {
      first = false;
    }
    prev_score = score;
  }
}


/***
 * \summary use the expectation maximization algorithm to estimate the
 *          parameters for this model given a set of sequences, their structure
 *          and diagnostic events.
 * \param seqs                TODO
 * \param secStr              TODO
 * \param diagEvents          diagnostic_events[i][j] is the location of the jth
 *                            diagnostic event in the ith sequence.
 *                            diagnostic_events[i] can be empty for any value
 *                            of i (i.e. there are no diagnostic events in that
 *                            sequence).
 * \param siteInd             TODO
 * \param seqInd              TODO
 */
void
Model::expectationMax_SeqStrDE(const vector<string> &seqs,
                               const vector<vector<double> > &secStr,
                               const vector<vector<double> > &diagEvents,
                               vector<vector<double> > &siteInd,
                               vector<double> &seqInd) {
  estimateDelta(seqs, diagEvents);
  double prev_score = std::numeric_limits<double>::max();
  double score = 0.0;
  bool first = true;
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_str_de(seqs, secStr, diagEvents, matrix, motif_sec_str,
                           f, f_sec_str, p, delta, gamma, siteInd, seqInd);
    maximization_seq(seqs, siteInd, seqInd, matrix, f, gamma);
    maximization_str(secStr, siteInd, seqInd, matrix, motif_sec_str, f_sec_str);
    maximization_de(diagEvents, siteInd, seqInd, matrix, p, delta);
    score = calculate_zoops_log_l(seqs, secStr, diagEvents, siteInd, seqInd);
    if (!first) {
      const double delta = fabs(prev_score - score);
      const double deltaProp = delta / fabs(prev_score);
      if (Model::DEBUG_LEVEL >= 1)
        cerr << "new ll: " << score << " old ll: " << prev_score << " "
             << "absolute change is " << delta << " fraction change "
             << deltaProp << endl;
      if (deltaProp < tolerance) break;
    } else {
      first = false;
    }
    prev_score = score;
  }
}


/***
 * \summary Estimate the parameters for this model using the expectation
 *          maximization algorithm. This function is a dispatcher; it selects
 *          the correct methods to call depending on whether the secondary
 *          structure of the sequences is known or not.
 * \param seqs          TODO
 * \param diagEvents    diagnostic_events[i][j] is the location of the jth
 *                      diagnostic event in the ith sequence.
 *                      diagnostic_events[i] can be empty for any value
 *                      of i (i.e. there are no diagnostic events in that
 *                      sequence).
 * \param secStruct     This vector must be either empty (if the secondary
 *                      structure is not known), or secStruct[i][j] gives the
 *                      probability that the jth nuc. in the ith sequence is
 *                      double-stranded (paired).
 * \param site_indic    TODO
 * \param seq_indic     TODO
 */
void
Model::expectationMax(const vector<string> &seqs,
                      const vector<vector<double> > &diagEvents,
                      const vector<vector<double> > &secStruct,
                      vector<vector<double> > &site_indic,
                      vector<double> &seq_indic) {
  if (secStruct.empty()) {
    expectation_maximization_seq_de(seqs, diagEvents, site_indic, seq_indic);
  }
  else {
    // TODO -- throw an exception here if the size of the secStrc vectors
    // don't match the sequences...?
    expectationMax_SeqStrDE(seqs, secStruct, diagEvents, site_indic, seq_indic);
  }
}


/******************************************************************************
 * SETTERS -- PARAMETER ESTIMATION BY EXHAUSTIVE SEARCH
 ******************************************************************************/

/***
 * \summary estimate the delta parameter for this model given a set of
 *          sequences and the diagnostic events within them.
 * \param seqs          seqs[i] is the ith sequence (string)
 * \param diagEvents    diagEvents[i][j] is the jth diagnostic event in the
 *                      ith sequence
 */
void
Model::estimateDelta(const vector<string> &seqs,
                     const vector<vector<double> > &diagEvents) {
  // if we have no diagnostic events, then just set delta to be the default,
  // otherwise we just try an exhaustive search of reasonable options with
  // a uniform PWM.
  size_t numDiagEvents = 0;
  for (size_t i=0; i < diagEvents.size(); i++)
    numDiagEvents += diagEvents[i].size();
  if (numDiagEvents == 0) {
    this->delta = Model::DEFAULT_DELTA;
  } else {
    vector<double> ll_delta;
    for (int deltaP = Model::MIN_DELTA; deltaP <= Model::MAX_DELTA; ++deltaP) {
      Model dummy;
      Model::set_model_uniform(this->size(), dummy);
      dummy.delta = deltaP;

      vector<double> has_motif(seqs.size(), dummy.gamma);
      vector<vector<double> > indicators;
      for (size_t i = 0; i < seqs.size(); ++i) {
        const size_t n_pos = seqs[i].length() - dummy.matrix.size() + 1;
        indicators.push_back(vector<double>(n_pos, 1.0 / n_pos));
      }
      // use EM to fit the sequence component and the geometric distribution,
      // but hold delta fixed (to avoid cyclic dependency)
      dummy.expectation_maximization_seq_de(seqs, diagEvents, indicators,
                                            has_motif, Model::HOLD_DELTA_FIXED);
      double score = calculate_zoops_log_l(seqs, diagEvents, indicators, has_motif);
      ll_delta.push_back(score);
    }

    double max_ll = -1.0 * std::numeric_limits<double>::max();
    int max_i = 0;
    for (size_t i = 0; i < ll_delta.size(); i++) {
      if (ll_delta[i] > max_ll) {
        max_ll = ll_delta[i];
        max_i = i + Model::MIN_DELTA;
      }
    }

    this->delta = max_i;
  }
}


/******************************************************************************
 * STATIC CLASS MEMBER FUNCTIONS
 ******************************************************************************/

/***
 * \summary TODO
 */
void
Model::set_model_uniform(const size_t width,
                         Model &model) {
  model.matrix.clear();
  model.matrix.resize(width,
                      vector<double>(alphabet_size, 1.0 / alphabet_size));
  model.motif_sec_str = vector<double>(width, 0.5);
  model.f = vector<double>(alphabet_size, 1.0 / alphabet_size);
  model.f_sec_str = 0.5;
  model.p = 0.5;
  model.delta = 0;
  model.gamma = 0.5;
}

/***
 * \summary TODO
 */
void
Model::set_model_by_word(const double pseudocount,
                         const string &kmer,
                         Model &model) {

  // initialize the matrix
  const size_t len = kmer.length();
  model.matrix.clear();
  model.matrix.resize(len, vector<double>(alphabet_size, pseudocount));

  // set the matrix to the word
  for (size_t i = 0; i < len; ++i) {
    const size_t base = base2int(kmer[i]);
    assert(base < alphabet_size);
    model.matrix[i][base] += 1.0;
  }

  // normalize matrix columns
  for (size_t i = 0; i < len; ++i)
    for (size_t j = 0; j < alphabet_size; ++j) {
      const double tot = accumulate(model.matrix[i].begin(),
                                    model.matrix[i].end(), 0.0);
      transform(model.matrix[i].begin(), model.matrix[i].end(),
                model.matrix[i].begin(),
                std::bind2nd(std::divides<double>(), tot));
    }
}

