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

using smithlab::alphabet_size;

// initialization of non-integral constants
const double Model::pseudocount = 0.1;
const double Model::tolerance = 1e-10;
const double Model::zoops_threshold = 0;
const double Model::DEFAULT_GEO_P = 0.5;

/******************************************************************************
 *        STATIC HELPER FUNCTIONS USED FOR CHECKING DATA CONSISTENCY
 *****************************************************************************/

/***
 * \summary TODO
 */
static void
checkAndThrow_consistent(const vector<string> &seqs,
                         const vector<vector<double> > &struc,
                         const vector<vector<size_t> > &diagEvents,
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
 * \summary TODO
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
 * \summary TODO
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
 * \summary TODO
 */
static void
calculate_number_of_bases_fg_bg(const vector<string> &sequences,
                                const vector<vector<double> > &site_indic,
                                const size_t motif_width,
                                vector<vector<double> > &nb_fg,
                                vector<double> &nb_bg) {

  nb_fg.clear();
  nb_fg.resize(motif_width,
               vector<double>(alphabet_size, 0.0));
  for (size_t i = 0; i < site_indic.size(); ++i)
    for (size_t j = 0; j < site_indic[i].size(); ++j)
      for (size_t k = 0; k < motif_width; ++k)
        nb_fg[k][base2int(sequences[i][j + k])] += site_indic[i][j];

  nb_bg.clear();
  nb_bg.resize(alphabet_size, 0.0);
  for (size_t i = 0; i < sequences.size(); ++i)
    for (size_t j = 0; j < sequences[i].length(); ++j)
      nb_bg[base2int(sequences[i][j])] += motif_width;

  for (size_t i = 0; i < site_indic.size(); ++i)
    for (size_t j = 0; j < site_indic[i].size(); ++j)
      for (size_t k = 0; k < motif_width; ++k)
        nb_bg[base2int(sequences[i][j + k])] -= site_indic[i][j];
}

/***
 * \summary TODO
 */
static void
calculate_number_of_bases_fg_bg_str(const vector<string> &sequences,
                                    const vector<vector<double> > &secondary_structure,
                                    const vector<vector<double> > &site_indic,
                                    const size_t motif_width,
                                    vector<vector<double> > &nb_fg_ss,
                                    vector<double> &nb_bg_ss,
                                    vector<vector<double> > &nb_fg_ds,
                                    vector<double> &nb_bg_ds) {

  nb_fg_ss.clear();
  nb_fg_ss.resize(motif_width, vector<double>(alphabet_size, 0.0));
  nb_fg_ds.clear();
  nb_fg_ds.resize(motif_width, vector<double>(alphabet_size, 0.0));

  // this loop calculates the number of bases occurring at a
  // particular location of the motif with particular secondary
  // structure
  for (size_t i = 0; i < site_indic.size(); ++i)
    for (size_t j = 0; j < site_indic[i].size(); ++j)
      for (size_t k = 0; k < motif_width; ++k) {
        const char base_idx = base2int(sequences[i][j + k]);
        const double curr_sec_str = secondary_structure[i][j + k];
        const double curr_site_indic = site_indic[i][j];
        nb_fg_ss[k][base_idx] += (1.0 - curr_sec_str) * curr_site_indic;
        nb_fg_ds[k][base_idx] += curr_sec_str * curr_site_indic;
      }

  nb_bg_ss.clear();
  nb_bg_ss.resize(alphabet_size, 0.0);
  nb_bg_ds.clear();
  nb_bg_ds.resize(alphabet_size, 0.0);

  // the next two loop sets calculate the number of bases occurring in
  // the background with a particular secondary structure
  for (size_t i = 0; i < sequences.size(); ++i)
    for (size_t j = 0; j < sequences[i].length(); ++j) {
      nb_bg_ss[base2int(sequences[i][j])] += (1.0 - secondary_structure[i][j])
          * motif_width;
      nb_bg_ds[base2int(sequences[i][j])] += (secondary_structure[i][j]
          * motif_width);
    }
  for (size_t i = 0; i < site_indic.size(); ++i)
    for (size_t j = 0; j < site_indic[i].size(); ++j)
      for (size_t k = 0; k < motif_width; ++k) {
        nb_bg_ss[base2int(sequences[i][j + k])] -= (1.0
            - secondary_structure[i][j + k])
                                                   * site_indic[i][j];
        nb_bg_ds[base2int(sequences[i][j + k])] -=
            (secondary_structure[i][j + k] * site_indic[i][j]);
      }
}


/******************************************************************************
 *          STATIC HELPER FUNCTIONS USED DURING MODEL MAXIMIZATION
 ******************************************************************************/

/***
 * \summary TODO
 */
static void
maximization_seq(const vector<string> &sequences,
                 const vector<vector<double> > &site_indic,
                 vector<double> &seq_indic,
                 vector<vector<double> > &matrix,
                 vector<double> &freq,
                 double &gamma) {

  static const double pseudocount = 1e-6;
  static const double TINY = 1e-100;

  vector<vector<double> > nb_fg(matrix.size(),
  // This value needs to be changed
                                vector<double>(alphabet_size, pseudocount));
  vector<double> nb_bg(alphabet_size, pseudocount);
  calculate_number_of_bases_fg_bg(sequences, site_indic, matrix.size(), nb_fg,
                                  nb_bg);

  for (size_t i = 0; i < matrix.size(); ++i) {
    const double total = accumulate(nb_fg[i].begin(), nb_fg[i].end(), 0.0);
    transform(nb_fg[i].begin(), nb_fg[i].end(), matrix[i].begin(),
              std::bind2nd(std::divides<double>(), total));

    // don't let anything get to zero; that'll be bad for log later.
    for (size_t j = 0; j < matrix[i].size(); ++j)
      if (matrix[i][j] < TINY) matrix[i][j] = TINY;
  }

  const double total = accumulate(nb_bg.begin(), nb_bg.end(), 0.0);
  transform(nb_bg.begin(), nb_bg.end(), freq.begin(),
            std::bind2nd(std::divides<double>(), total));

  gamma = accumulate(seq_indic.begin(), seq_indic.end(), 0.0)
      / sequences.size();
}

/***
 * \summary maximize parameters for the structure component of the model based
 *          on the data and the indicators.
 * \param sequences             TODO
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
maximization_str(const vector<string> &sequences,
                 const vector<vector<double> > &secondary_structure,
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

/***
 * \summary TODO
 */
static void
maximization_de(const vector<string> &sequences,
                const vector<vector<size_t> > &diagEvents,
                const vector<vector<double> > &site_indic,
                const vector<double> &seq_indic,
                vector<vector<double> > &matrix,
                double &geo_p,
                int &geo_delta) {
  // if we have no diagnostic events, we can't do anything.. just set this to
  // some default value
  size_t numDEs = 0;
  for (size_t i = 0; i < diagEvents.size(); ++i) numDEs += diagEvents[i].size();
  if (numDEs == 0) {
    geo_p = Model::DEFAULT_GEO_P;
  } else {
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < site_indic.size(); i++) {
      double seq_sum = 0.0;
      for (size_t k = 0; k < site_indic[i].size(); k++) {
        numerator += (seq_indic[i] * site_indic[i][k] * diagEvents[i].size());
        double site_sum = 0.0;
        if (diagEvents[i].size() > 0)
          for (size_t j = 0; j < diagEvents[i].size(); j++)
            site_sum += abs(diagEvents[i][j] - (k + geo_delta));
        seq_sum += site_indic[i][k] * (diagEvents[i].size() + site_sum);
      }
      denominator += seq_indic[i] * seq_sum;
    }
    geo_p = max(min(numerator / denominator, 0.999),
                std::numeric_limits<double>::min());
  }
}


/******************************************************************************
 *   STATIC HELPER FUNCTIONS USED IN CALCULATION OF INDICATORS / EXPECTATION
 ******************************************************************************/

/***
 * \summary TODO
 */
static void
get_numerator_seq_de_for_site(const string &seq,
                              const vector<size_t> &diagnostic_events,
                              const vector<vector<double> > &matrix,
                              const vector<double> &freqs,
                              const double geo_p,
                              const int geo_delta,
                              const double gamma,
                              const size_t site,
                              double &num) {

  vector<double> f_powers(alphabet_size, 0.0);
  for (size_t i = 0; i < seq.length(); ++i) {
    const size_t base = base2int(seq[i]);
    if (i >= site && i < site + matrix.size())
      num += log(matrix[i - site][base]);
    else
      f_powers[base]++;
    if (!std::isfinite(f_powers[base]))
      throw SMITHLABException("failed expectation calc; f_powers non-finite");
    if (!std::isfinite(num)) {
      stringstream ss;
      ss << "failed expectation calculation; numerator non-finite. Matrix "
         << "entry was: " << matrix[i - site][base];
      throw SMITHLABException(ss.str());
    }
  }
  for (size_t b = 0; b < alphabet_size; b++)
    num += f_powers[b] * log(freqs[b]);
  if (diagnostic_events.size() > 0) {
    double power = 0.0;
    for (size_t j = 0; j < diagnostic_events.size(); j++)
      power += abs(diagnostic_events[j] - (site + geo_delta));
    num += ((power * log(1 - geo_p)) + (diagnostic_events.size() * log(geo_p)));
  }

  num += log(gamma / (seq.length() - matrix.size() + 1.0));
}

/***
 * \summary TODO
 */
static void
get_numerator_seq_str_de_for_site(const string &seq,
                               const vector<double> &secondary_structure,
                               const vector<size_t> &diagnostic_events,
                               const vector<vector<double> > &matrix,
                               const vector<double> &motif_sec_str,
                               const vector<double> &freqs,
                               const double f_sec_str,
                               const double geo_p,
                               const int geo_delta,
                               const double gamma,
                               const size_t site,
                               double &num) {

  // calculating the contribution of the foreground and the powers
  // that will be needed for the background calculations (below).
  vector<double> f_powers_ss(alphabet_size, 0.0);
  vector<double> f_powers_ds(alphabet_size, 0.0);
  for (size_t i = 0; i < seq.length(); ++i) {
    const size_t base = base2int(seq[i]);
    if (i >= site && i < site + matrix.size()) {
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
    } else {
      f_powers_ss[base] += (1.0 - secondary_structure[i]);
      f_powers_ds[base] += (secondary_structure[i]);
    }
    if (!std::isfinite(f_powers_ss[base]))
      throw SMITHLABException("failed expectation calc; f_powers_ss non-finite");
    if (!std::isfinite(f_powers_ds[base]))
      throw SMITHLABException("failed expectation calc; f_powers_ds non-finite");
  }

  // calculating the contribution of the background (outside the motif
  // occurrences)
  for (size_t b = 0; b < alphabet_size; b++) {
    num += f_powers_ss[b] * log(freqs[b] * (1 - f_sec_str));
    num += f_powers_ds[b] * log(freqs[b] * (f_sec_str));
  }
  if (diagnostic_events.size() > 0) {
    double power = 0.0;
    for (size_t j = 0; j < diagnostic_events.size(); j++)
      power += abs(diagnostic_events[j] - (site + geo_delta));
    num += ((power * log(1 - geo_p)) + (diagnostic_events.size() * log(geo_p)));
  }
  num += log(gamma / (seq.length() - matrix.size() + 1.0));
}

/***
 * \summary TODO
 */
static void
expectation_seq_str_de_for_single_seq(const string &seq,
                                      const vector<double> &secondary_structure,
                                      const vector<size_t> &diagnostic_events,
                                      const vector<vector<double> > &matrix,
                                      const vector<double> &motif_sec_str,
                                      const vector<double> &freqs,
                                      const double f_sec_str,
                                      const double geo_p,
                                      const int geo_delta,
                                      const double gamma,
                                      vector<double> &site_indic,
                                      double &seq_indic) {


  // get log likelihood for each site
  vector<double> numerator(site_indic.size(), 0.0);
  for (size_t i = 0; i < site_indic.size(); ++i)
    get_numerator_seq_str_de_for_site(seq, secondary_structure,
                                      diagnostic_events, matrix, motif_sec_str,
                                      freqs, f_sec_str, geo_p, geo_delta, gamma,
                                      i, numerator[i]);

  double no_motif = 0.0;
  for (size_t i = 0; i < seq.length(); i++) {
    no_motif += secondary_structure[i]
        * log(freqs[base2int(seq[i])] * f_sec_str);
    no_motif += ((1 - secondary_structure[i])
        * log(freqs[base2int(seq[i])] * (1 - f_sec_str)));
    if (diagnostic_events.size() > 0)
      no_motif -= (diagnostic_events.size() * log(seq.length()));
  }
  numerator.push_back(no_motif + log(1.0 - gamma));

  const double denominator = smithlab::log_sum_log_vec(numerator,
                                                       numerator.size());
  for (size_t i = 0; i < site_indic.size(); ++i)
    site_indic[i] = exp(numerator[i] - denominator);

  seq_indic = accumulate(site_indic.begin(), site_indic.end(), 0.0);
}

/***
 * \summary TODO
 */
static void
expectation_seq_de_for_single_seq(const string &seq,
                                  const vector<size_t> &diagnostic_events,
                                  const vector<vector<double> > &matrix,
                                  const vector<double> &freqs,
                                  const double geo_p,
                                  const int geo_delta,
                                  const double gamma,
                                  vector<double> &site_indic,
                                  double &seq_indic) {
  // get log likelihood for each site
  vector<double> numerator(site_indic.size(), 0.0);
  for (size_t i = 0; i < site_indic.size(); ++i)
    get_numerator_seq_de_for_site(seq, diagnostic_events, matrix, freqs, geo_p,
                                  geo_delta, gamma, i, numerator[i]);

  double no_motif = 0.0;
  for (size_t i = 0; i < seq.length(); i++)
    no_motif += log(freqs[base2int(seq[i])]);
  if (diagnostic_events.size() > 0)
    no_motif -= (diagnostic_events.size() * log(seq.length()));
  numerator.push_back(no_motif + log(1.0 - gamma));

  const double denominator = smithlab::log_sum_log_vec(numerator,
                                                       numerator.size());
  for (size_t i = 0; i < site_indic.size(); ++i)
    site_indic[i] = exp(numerator[i] - denominator);

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
  vector<size_t> diagEvents;
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
                   const vector<vector<size_t> > &diagnostic_events,
                   const vector<vector<double> > &matrix,
                   const vector<double> &freqs,
                   const double geo_p,
                   const int geo_delta,
                   const double gamma,
                   vector<vector<double> > &site_indic,
                   vector<double> &seq_indic) {
  if (Model::DEBUG_MESSAGES)
    cerr << "performing expectation step with matrix " << endl
         << matrixToString(matrix) << endl;
  for (size_t i = 0; i < sequences.size(); i++)
    expectation_seq_de_for_single_seq(sequences[i], diagnostic_events[i],
                                      matrix, freqs, geo_p, geo_delta, gamma,
                                      site_indic[i], seq_indic[i]);
  if (Model::DEBUG_MESSAGES)
    cerr << "finished expectation step" << endl;
}

/***
 * \summary TODO
 */
static void
expectation_seq_str_de(const vector<string> &sequences,
                       const vector<vector<double> > &secondary_structure,
                       const vector<vector<size_t> > &diagnostic_events,
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
 * GETTERS
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
  vector<vector<size_t> > diagEvents (seqs.size());
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
                             const vector<vector<size_t> > &diagnostic_events,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {

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

  for (size_t i = 0; i < sequences.size(); i++) {
    if (diagnostic_events[i].size() > 0) {
      for (size_t k = 0; k < site_indic[i].size(); k++) {
        double power = 0.0;
        for (size_t j = 0; j < diagnostic_events[i].size(); j++) {
          power += abs(diagnostic_events[i][j] - (k + delta));
          if (!std::isfinite(power)) {
            stringstream ss;
            ss << "failed likelihood calculation; power is non-finite. Diag. "
               << "event is: " << diagnostic_events[i][j] << "; k+delta is: "
               << (k+delta);
            throw SMITHLABException(ss.str());
          }
        }
        ret += site_indic[i][k]
            * ((power * log(1 - p)) + (diagnostic_events[i].size() * log(p)));
      }
    }
  }

  //--------
  for (size_t i = 0; i < sequences.size(); i++) {
    double has_no_motif = 0.0;
    for (size_t j = 0; j < sequences[i].length(); j++)
      has_no_motif += log(f[base2int(sequences[i][j])]);
    if (diagnostic_events[i].size() > 0)
      has_no_motif -= (diagnostic_events[i].size() * log(sequences[i].length()));
    ret += (1 - seq_indic[i]) * has_no_motif;
    ret += (1 - seq_indic[i]) * log(1 - gamma);
    ret += seq_indic[i]
        * log(gamma / site_indic[i].size());
  }
  //--------

  return ret;
}

/***
 * \summary calculate the log-likelihood for this model given a set of
 *          sequences, their structure, and indicator variables on occurrences
 *          of the motif
 * \param sequences             TODO
 * \param secondary_structure   secondary_structure[i][j] is the probability
 *                              that the jth base in the ith sequence is paired.
 * \param site_indic            site_indic[i][j] is the probability that the
 *                              jth position in the ith sequence is the start
 *                              of a motif occurrence.
 * \param seq_indic             seq_indic[i] is the probability that the ith
 *                              sequence contains an occurrence of the motif.
 */
double
Model::calculate_zoops_log_l(const vector<string> &sequences,
                             const vector<vector<double> > &secondary_structure,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {
  // calc. without diag. events is equivalent to having no diag. events, so...
  vector< vector<size_t> > diagEvents(sequences.size());
  return calculate_zoops_log_l(sequences, secondary_structure, diagEvents,
                        site_indic, seq_indic);
}

/***
 * \summary calculate the zoops log-likelihood for this model given a set of
 *          sequences, secondary structure and diagnostic events.
 */
double
Model::calculate_zoops_log_l(const vector<string> &sequences,
                             const vector<vector<double> > &secondary_structure,
                             const vector<vector<size_t> > &diagnostic_events,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {
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
    ret += (nb_bg_ss[i] * log(f[i] * (1 - f_sec_str)));
    ret += (nb_bg_ds[i] * log(f[i] * f_sec_str));
    for (size_t j = 0; j < matrix.size(); ++j) {
      ret += (nb_fg_ss[j][i] * log(matrix[j][i] * (1 - motif_sec_str[j])));
      ret += (nb_fg_ds[j][i] * log(matrix[j][i] * motif_sec_str[j]));
    }
  }

  // ZOOPS-specific calculation starts from here --------
  for (size_t i = 0; i < sequences.size(); i++) {
    if (diagnostic_events[i].size() > 0) {
      for (size_t k = 0; k < site_indic[i].size(); k++) {
        double power = 0.0;
        for (size_t j = 0; j < diagnostic_events[i].size(); j++) {
          power += abs(diagnostic_events[i][j] - (k + delta));
          if (!std::isfinite(power)) {
            stringstream ss;
            ss << "failed likelihood calculation; power is non-finite. Diag. "
               << "event is: " << diagnostic_events[i][j] << "; k+delta is: "
               << (k+delta);
            throw SMITHLABException(ss.str());
          }
        }
        ret += site_indic[i][k]
            * ((power * log(1 - p)) + (diagnostic_events[i].size() * log(p)));
      }
    }
  }

  for (size_t i = 0; i < sequences.size(); i++) {
    double has_no_motif = 0.0;
    for (size_t j = 0; j < sequences[i].length(); j++)
      has_no_motif += log(f[base2int(sequences[i][j])]);
    if (diagnostic_events[i].size() > 0)
      has_no_motif -=
          (diagnostic_events[i].size() * log(sequences[i].length()));
    ret += (1 - seq_indic[i]) * has_no_motif;
    ret += (1 - seq_indic[i]) * log(1 - gamma);
    ret += seq_indic[i] * log(gamma / site_indic[i].size());
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
  vector<vector<size_t> > diagEvents(sequences.size());
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
      if (Model::DEBUG_MESSAGES)
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
                                       const vector<vector<size_t> > &diagEvents,
                                       vector<vector<double> > &site_indic,
                                       vector<double> &seq_indic,
                                       const bool holdDelta = false) {
  if (!holdDelta) this->estimateDelta(seqs, diagEvents);
  double prev_score = std::numeric_limits<double>::max();
  bool first = true;
  double score = 0.0;
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_de(seqs, diagEvents, matrix, f, p, delta, gamma,
                       site_indic, seq_indic);
    maximization_seq(seqs, site_indic, seq_indic, matrix, f, gamma);
    maximization_de(seqs, diagEvents, site_indic, seq_indic, matrix, p, delta);
    score = calculate_zoops_log_l(seqs, diagEvents, site_indic, seq_indic);
    if (!first) {
      const double delta = fabs(prev_score - score);
      const double deltaProp = delta / fabs(prev_score);
      if (Model::DEBUG_MESSAGES)
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
 * \param sequences           TODO
 * \param secStructure        TODO
 * \param diagnostic_events   diagnostic_events[i][j] is the location of the jth
 *                            diagnostic event in the ith sequence.
 *                            diagnostic_events[i] can be empty for any value
 *                            of i (i.e. there are no diagnostic events in that
 *                            sequence).
 * \param site_indic          TODO
 * \param seq_indic           TODO
 */
void
Model::expectationMax_SeqStrDE(const vector<string> &sequences,
                               const vector<vector<double> > &secStructure,
                               const vector<vector<size_t> > &diagnostic_events,
                               vector<vector<double> > &site_indic,
                               vector<double> &seq_indic) {
  estimateDelta(sequences, diagnostic_events);
  double prev_score = std::numeric_limits<double>::max();
  double score = 0.0;
  bool first = true;
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_str_de(sequences, secStructure, diagnostic_events,
                           matrix, motif_sec_str, f, f_sec_str, p, delta, gamma,
                           site_indic, seq_indic);
    maximization_seq(sequences, site_indic, seq_indic, matrix, f, gamma);
    maximization_str(sequences, secStructure, site_indic, seq_indic,
                     matrix, motif_sec_str, f_sec_str);
    maximization_de(sequences, diagnostic_events, site_indic, seq_indic, matrix,
                    p, delta);
    score = calculate_zoops_log_l(sequences, secStructure,
                                  diagnostic_events, site_indic, seq_indic);
    if (!first) {
      const double delta = fabs(prev_score - score);
      const double deltaProp = delta / fabs(prev_score);
      if (Model::DEBUG_MESSAGES)
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
                      const vector<vector<size_t> > &diagEvents,
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
                     const vector<vector<size_t> > &diagEvents) {
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
        max_i = i - 10;
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
  for (size_t i = 0; i < len; ++i)
    model.matrix[i][base2int(kmer[i])] += 1.0;

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

