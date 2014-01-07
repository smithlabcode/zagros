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

// TO HELP WITH DEBUGGING:
// static string
// format_matrix(const vector<vector<double> > &matrix) {
//   std::ostringstream oss;
//   for (size_t i = 0; i < matrix.size(); ++i) {
//     copy(matrix[i].begin(), matrix[i].end(), 
// 	 std::ostream_iterator<double>(oss, "\t"));
//     oss << endl;
//   }
//   oss << endl;
//   return oss.str();
// }

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

void
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

double
Model::calculate_zoops_log_l(const vector<string> &sequences,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {

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

  //--------
  for (size_t i = 0; i < sequences.size(); i++) {
    double has_no_motif = 0.0;
    for (size_t j = 0; j < sequences[i].length(); j++)
      has_no_motif += log(f[base2int(sequences[i][j])]);
    ret += (1 - seq_indic[i]) * has_no_motif;
    ret += (1 - seq_indic[i]) * log(1 - gamma);
    ret += seq_indic[i]
        * log(gamma / site_indic.size());
  }
  //--------

  return ret;
}

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

void
Model::expectation_maximization(const vector<string> &sequences,
                                const vector<vector<size_t> > &diagnostic_events,
                                const vector<vector<double> > &secondary_structure,
                                vector<vector<double> > &site_indic,
                                vector<double> &seq_indic) {
  if (secondary_structure.empty() && diagnostic_events.empty())
    expectation_maximization_seq(sequences, site_indic, seq_indic);
  else if (secondary_structure.empty())
    expectation_maximization_seq_de(sequences, diagnostic_events, site_indic,
                                    seq_indic);
  else if (diagnostic_events.empty())
    expectation_maximization_seq_str(sequences, secondary_structure, site_indic,
                                     seq_indic);
  else
    expectation_maximization_seq_str_de(sequences, secondary_structure,
                                        diagnostic_events, site_indic, seq_indic);
}


static void
get_numerator_for_site(const string &seq,
                       const vector<vector<double> > &matrix,
                       const vector<double> &freqs,
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

    /*using std::cout;
    using std::endl;
    cout << "in get num" << endl;
    for (size_t i = 0; i < matrix.size(); ++i) {
      for (size_t j = 0; j < matrix[i].size(); ++j) {
        cout << matrix[i][j] << ", ";
      }
      cout << endl;
    }*/

    assert(std::isfinite(f_powers[base]) && std::isfinite(num));
  }
  for (size_t b = 0; b < alphabet_size; b++)
    num += f_powers[b] * log(freqs[b]);
  num += log(gamma / (seq.length() - matrix.size() + 1.0));
}

static void
expectation_for_single_seq(const string &seq,
                           const vector<vector<double> > &matrix,
                           const vector<double> &freqs,
                           const double gamma,
                           vector<double> &site_indic,
                           double &seq_indic) {

  // get log likelihood for each site
  vector<double> numerator(site_indic.size(), 0.0);
  for (size_t i = 0; i < site_indic.size(); ++i)
    get_numerator_for_site(seq, matrix, freqs, gamma, i, numerator[i]);

  double no_motif = 0.0;
  for (size_t i = 0; i < seq.length(); i++)
    no_motif += log(freqs[base2int(seq[i])]);
  numerator.push_back(no_motif + log(1.0 - gamma));

  const double denominator = smithlab::log_sum_log_vec(
      numerator, numerator.size());
  for (size_t i = 0; i < site_indic.size(); ++i)
    site_indic[i] = exp(numerator[i] - denominator);

  seq_indic = accumulate(site_indic.begin(), site_indic.end(), 0.0);
}

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
    for (size_t j = 0; j < matrix[i].size(); ++j) if (matrix[i][j] < TINY) matrix[i][j] = TINY;
  }

  const double total = accumulate(nb_bg.begin(), nb_bg.end(), 0.0);
  transform(nb_bg.begin(), nb_bg.end(), freq.begin(),
            std::bind2nd(std::divides<double>(), total));

  gamma = accumulate(seq_indic.begin(), seq_indic.end(), 0.0)
      / sequences.size();

  using std::cout;
  using std::endl;
  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix[i].size(); ++j) {
      cout << matrix[i][j] << ", ";
    }
    cout << endl;
  }
}

void
Model::expectation_maximization_seq(const vector<string> &sequences,
                                    vector<vector<double> > &site_indic,
                                    vector<double> &seq_indic) {

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {

    expectation_seq(sequences, matrix, f, gamma, site_indic, seq_indic);
    maximization_seq(sequences, site_indic, seq_indic, matrix, f, gamma);

    const double score = calculate_zoops_log_l(sequences, site_indic,
                                               seq_indic);

    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////
///////////  CODE FOR SEQUENCE AND STRUCTURE
///////////

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

double
Model::calculate_zoops_log_l(const vector<string> &sequences,
                             const vector<vector<double> > &secondary_structure,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {

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

  // this is the ZOOPS-specific part of the log likelihood
  for (size_t i = 0; i < sequences.size(); i++) {
    double has_no_motif = 0.0;
    for (size_t j = 0; j < sequences[i].length(); j++)
      has_no_motif += log(f[base2int(sequences[i][j])]);
    ret += (1.0 - seq_indic[i]) * has_no_motif;
    ret += (1.0 - seq_indic[i]) * log(1.0 - gamma);
    ret += seq_indic[i]
        * log(gamma / (sequences[i].length() - matrix.size() + 1));
  }

  return ret;
}

static void
get_numerator_seq_str_for_site(const string &seq,
                               const vector<double> &secondary_structure,
                               const vector<vector<double> > &matrix,
                               const vector<double> &motif_sec_str,
                               const vector<double> &freqs,
                               const double f_sec_str,
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
    } else {
      f_powers_ss[base] += (1.0 - secondary_structure[i]);
      f_powers_ds[base] += (secondary_structure[i]);
    }
    assert(std::isfinite(f_powers_ss[base]) && std::isfinite(f_powers_ds[base])
           && std::isfinite(num));
  }

  // calculating the contribution of the background (outside the motif
  // occurrences)
  for (size_t b = 0; b < alphabet_size; b++) {
    num += f_powers_ss[b] * log(freqs[b] * (1 - f_sec_str));
    num += f_powers_ds[b] * log(freqs[b] * (f_sec_str));
  }
  num += log(gamma / (seq.length() - matrix.size() + 1.0));
}

static void
expectation_seq_str_for_single_seq(const string &seq,
                                   const vector<double> &secondary_structure,
                                   const vector<vector<double> > &matrix,
                                   const vector<double> &motif_sec_str,
                                   const vector<double> &freqs,
                                   const double f_sec_str,
                                   const double gamma,
                                   vector<double> &site_indic,
                                   double &seq_indic) {

  // get log likelihood for each site
  vector<double> numerator(site_indic.size(), 0.0);
  for (size_t i = 0; i < site_indic.size(); ++i)
    get_numerator_seq_str_for_site(seq, secondary_structure, matrix,
                                   motif_sec_str, freqs, f_sec_str, gamma, i,
                                   numerator[i]);

  double no_motif = 0.0;
  for (size_t i = 0; i < seq.length(); i++) {
    no_motif += secondary_structure[i]
        * log(freqs[base2int(seq[i])] * f_sec_str);
    no_motif += ((1.0 - secondary_structure[i])
        * log(freqs[base2int(seq[i])] * (1.0 - f_sec_str)));
  }
  numerator.push_back(no_motif + log(1.0 - gamma));

  const double denominator = smithlab::log_sum_log_vec(numerator,
                                                       numerator.size());
  for (size_t i = 0; i < site_indic.size(); ++i)
    site_indic[i] = exp(numerator[i] - denominator);

  seq_indic = accumulate(site_indic.begin(), site_indic.end(), 0.0);
}

static void
expectation_seq_str(const vector<string> &sequences,
                    const vector<vector<double> > &secondary_structure,
                    const vector<vector<double> > &matrix,
                    const vector<double> &motif_sec_str,
                    const vector<double> &freqs,
                    const double f_sec_str,
                    const double gamma,
                    vector<vector<double> > &site_indic,
                    vector<double> &seq_indic) {

  for (size_t i = 0; i < sequences.size(); i++) {
      expectation_seq_str_for_single_seq(sequences[i], secondary_structure[i],
                                         matrix, motif_sec_str, freqs, f_sec_str,
                                         gamma, site_indic[i], seq_indic[i]);
  }
}

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
      for (size_t site = 0; site < site_indic[i].size(); ++site)
        motif_sec_str[i] += seq_indic[j] * site_indic[j][site]
                            * secondary_structure[j][site + i];
    }
    motif_sec_str[i] = motif_sec_str[i]
        / accumulate(seq_indic.begin(), seq_indic.end(), 0.0);
  }
  //If we want to learn this parameter we have to calculate this:
  f_sec_str = 0.5;
}

void
Model::expectation_maximization_seq_str(const vector<string> &sequences,
                                        const vector<vector<double> > &sec_structure,
                                        vector<vector<double> > &site_indic,
                                        vector<double> &seq_indic) {

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_str(sequences, sec_structure, matrix, motif_sec_str, f,
                        f_sec_str, gamma, site_indic, seq_indic);
    maximization_seq(sequences, site_indic, seq_indic, matrix, f, gamma);
    maximization_str(sequences, sec_structure, site_indic, seq_indic, matrix,
                     motif_sec_str, f_sec_str);

    const double score = calculate_zoops_log_l(sequences, sec_structure,
                                               site_indic, seq_indic);
    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////
///////////  CODE FOR SEQUENCE AND DIAGNOSTIC EVENTS
///////////

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
    assert(std::isfinite(f_powers[base]) && std::isfinite(num));
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
  for (size_t i = 0; i < sequences.size(); i++)
    expectation_seq_de_for_single_seq(sequences[i], diagnostic_events[i],
                                      matrix, freqs, geo_p, geo_delta, gamma,
                                      site_indic[i], seq_indic[i]);
}

static void
maximization_de(const vector<string> &sequences,
                const vector<vector<size_t> > &diagnostic_events,
                const vector<vector<double> > &site_indic,
                const vector<double> &seq_indic,
                vector<vector<double> > &matrix,
                double &geo_p,
                int &geo_delta) {

  double numerator = 0.0;
  double denominator = 0.0;
  for (size_t i = 0; i < site_indic.size(); i++) {
    double seq_sum = 0.0;
    for (size_t k = 0; k < site_indic[i].size(); k++) {
      numerator += (seq_indic[i] * site_indic[i][k] * diagnostic_events[i].size());
      double site_sum = 0.0;
      if (diagnostic_events[i].size() > 0)
        for (size_t j = 0; j < diagnostic_events[i].size(); j++)
          site_sum += abs(diagnostic_events[i][j] - (k + geo_delta));
      seq_sum += site_indic[i][k] * ( diagnostic_events[i].size() + site_sum);
    }
    denominator += seq_indic[i] * seq_sum;
  }
  geo_p = max(min(numerator / denominator, 0.999), std::numeric_limits<double>::min());
}

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
    for (size_t j = 0; j < matrix.size(); ++j)
      ret += nb_fg[j][i] * log(matrix[j][i]);
  }

  for (size_t i = 0; i < sequences.size(); i++) {
    if (diagnostic_events[i].size() > 0) {
      for (size_t k = 0; k < site_indic[i].size(); k++) {
        double power = 0.0;
        for (size_t j = 0; j < diagnostic_events[i].size(); j++)
          power += abs(diagnostic_events[i][j] - (k + delta));
        assert(std::isfinite(power));
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

static int
find_delta(const vector<string> &sequences,
           const vector<vector<size_t> > &diagnostic_events,
           const vector<vector<double> > &matrix) {

  vector<double> ll_delta;
  for (int delta_param = -10; delta_param <= 10; ++delta_param) {
    cerr << "\r" << int(100 * ll_delta.size() / 21) << "% completed..."
         << std::flush;

    Model m;
    Model::set_model_uniform(matrix.size(), m);
    m.delta = delta_param;

    vector<double> has_motif(sequences.size(), m.gamma);
    vector<vector<double> > indicators;
    for (size_t i = 0; i < sequences.size(); ++i) {
      const size_t n_pos = sequences[i].length() - m.matrix.size() + 1;
      indicators.push_back(vector<double>(n_pos, 1.0 / n_pos));
    }
    double prev_score = std::numeric_limits<double>::max();
    double score = 0.0;
    for (size_t i = 0; i < m.max_iterations; ++i) {
      expectation_seq_de(sequences, diagnostic_events, m.matrix, m.f, m.p,
                         m.delta, m.gamma, indicators, has_motif);
      maximization_seq(sequences, indicators, has_motif, m.matrix, m.f,
                       m.gamma);
      maximization_de(sequences, diagnostic_events, indicators, has_motif,
                      m.matrix, m.p, m.delta);

      score = m.calculate_zoops_log_l(sequences, diagnostic_events,
                                                   indicators, has_motif);
      if (abs(prev_score - score) / prev_score < Model::tolerance) {
        break;
      }
    }
    ll_delta.push_back(score);
  }
  cerr << "\r" << "100% completed..." << endl;

  double max_ll = -1.0 * std::numeric_limits<double>::max();
  int max_i = 0;
  for (size_t i = 0; i < ll_delta.size(); i++) {
    if (ll_delta[i] > max_ll) {
      max_ll = ll_delta[i];
      max_i = i - 10;
    }
  }
  return max_i;
}


void
Model::expectation_maximization_seq_de(const vector<string> &sequences,
                                       const vector<vector<size_t> > &diagnostic_events,
                                       vector<vector<double> > &site_indic,
                                       vector<double> &seq_indic) {

  delta = find_delta(sequences, diagnostic_events, matrix);
  double prev_score = std::numeric_limits<double>::max();
  double score = 0.0;
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_de(sequences, diagnostic_events, matrix, f, p, delta, gamma,
                       site_indic, seq_indic);
    maximization_seq(sequences, site_indic, seq_indic, matrix, f, gamma);
    maximization_de(sequences, diagnostic_events, site_indic, seq_indic, matrix,
                    p, delta);

    score = calculate_zoops_log_l(sequences, diagnostic_events,
                                               site_indic, seq_indic);
    if (abs(prev_score - score) / prev_score < tolerance) {
      break;
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////
///////////  CODE FOR SEQUENCE, STRUCTURE AND DIAGNOSTIC EVENTS
///////////


double
Model::calculate_zoops_log_l(const vector<string> &sequences,
                             const vector<vector<double> > &secondary_structure,
                             const vector<vector<size_t> > &diagnostic_events,
                             const vector<vector<double> > &site_indic,
                             const vector<double> &seq_indic) const {

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

  //--------
  for (size_t i = 0; i < sequences.size(); i++) {
    if (diagnostic_events[i].size() > 0) {
      for (size_t k = 0; k < site_indic[i].size(); k++) {
        double power = 0.0;
        for (size_t j = 0; j < diagnostic_events[i].size(); j++)
          power += abs(diagnostic_events[i][j] - (k + delta));
        assert(std::isfinite(power));
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
  //--------

  return ret;
}

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
    } else {
      f_powers_ss[base] += (1 - secondary_structure[i]);
      f_powers_ds[base] += (secondary_structure[i]);
    }
    assert(
        std::isfinite(f_powers_ss[base]) && std::isfinite(f_powers_ds[base])
        && std::isfinite(num));
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


void
Model::expectation_maximization_seq_str_de(const vector<string> &sequences,
                                           const vector<vector<double> > &secondary_structure,
                                           const vector<vector<size_t> > &diagnostic_events,
                                           vector<vector<double> > &site_indic,
                                           vector<double> &seq_indic) {
  delta = find_delta(sequences, diagnostic_events, matrix);
  double prev_score = std::numeric_limits<double>::max();
  double score = 0.0;
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_str_de(sequences, secondary_structure, diagnostic_events,
                           matrix, motif_sec_str, f, f_sec_str, p, delta, gamma,
                           site_indic, seq_indic);
    maximization_seq(sequences, site_indic, seq_indic, matrix, f, gamma);
    maximization_str(sequences, secondary_structure, site_indic, seq_indic,
                     matrix, motif_sec_str, f_sec_str);
    maximization_de(sequences, diagnostic_events, site_indic, seq_indic, matrix,
                    p, delta);

    score = calculate_zoops_log_l(sequences, secondary_structure,
                                  diagnostic_events, site_indic, seq_indic);
    if (abs(prev_score - score) / prev_score < tolerance) {
      break;
    }
  }
}
