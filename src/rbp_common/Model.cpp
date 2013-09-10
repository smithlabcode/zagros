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
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include <numeric>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <queue>
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Model.hpp"
#include "IO.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::accumulate;

using smithlab::alphabet_size;

Model::Model(const size_t motif_width) {

  for (size_t i = 0; i < motif_width; ++i)
    M.push_back(vector<double>(alphabet_size, 1.0 / alphabet_size));
  lambda.resize(motif_width, 0.5);
  f.resize(alphabet_size, 1.0 / alphabet_size);
  p = 0.5;
  delta = 0;
  gamma = 1;
}

double Model::calculateLogL(const vector<string> &sequences,
    const vector<vector<double> > &indicators, const vector<double> &zoops_i) {

  const size_t n = sequences.size();
  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0.0));

  for (size_t i = 0; i < n; i++) {
    const size_t l = sequences[i].length();
    for (size_t k = 0; k < indicators[i].size(); k++)
      for (size_t j = 0; j < l; j++)
        if (j >= k && j < k + M.size())
          nb[j - k + 1][base2int(sequences[i][j])] += indicators[i][k];
        else
          nb[0][base2int(sequences[i][j])] += indicators[i][k];
  }

  double ret = 0.0;
  for (size_t b = 0; b < alphabet_size; b++) {
    ret += nb[0][b] * log(f[b]);
    for (size_t j = 0; j < M.size(); j++)
      ret += nb[j + 1][b] * log(M[j][b]);
  }

  return ret;
}

void Model::expectation_maximization(const vector<string> &sequences,
    const vector<vector<size_t> > &diagnostic_events,
    const vector<double> &secondary_structure,
    vector<vector<double> > &indicators, vector<double> &zoops_i,
    const string t, const bool d, const string &file_name_base,
    const size_t max_iterations, const double tolerance) {

  vector<kmer_info> top_five_kmers;
//  determineStartingPoint_best_kmer(sequences, top_five_kmers);
  string starting_point = "GGCCCGCG";//top_five_kmers.front().kmer;
  cout << starting_point << endl;
  set_model(starting_point);
  expectation_maximization_seq(
      sequences, indicators, zoops_i, max_iterations, tolerance);
}

void Model::expectation_maximization_seq(const vector<string> &sequences,
    vector<vector<double> > &indicators, vector<double> &zoops_i,
    const size_t max_iterations, const double tolerance) {

  cerr << "Fitting the full model using sequence started...";

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq(sequences, indicators, zoops_i);
    maximization_seq(sequences, indicators, zoops_i);

    const double score = calculateLogL(sequences, indicators, zoops_i);

    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
    prev_score = score;
  }

  cerr << "done!" << endl;
}

void Model::expectation_seq(const vector<string> &sequences,
    vector<vector<double> > &indicators, vector<double> &zoops_i) {

  const size_t n = sequences.size();
  for (size_t i = 0; i < n; i++) {

    const size_t l = sequences[i].length();

    vector<double> numerator(indicators[i].size(), 0.0);
    //--------
    double has_motif = 0.0;
    //--------
    for (size_t k = 0; k < indicators[i].size(); k++) {

      vector<double> f_powers(alphabet_size, 0.0);
      for (size_t j = 0; j < l; j++) {
        const char base = base2int(sequences[i][j]);
        //--------
        has_motif += log(f[base]);
        //--------
        if (j >= k && j < k + M.size()) {
          numerator[k] += log(M[j - k][base]);
        } else
          f_powers[base2int(sequences[i][j])]++;
        assert(std::isfinite(f_powers[base]));
        if (!std::isfinite(numerator[k])) {
          cerr << "Numerical error!" << endl;
          exit(1);
        }
      }
      for (size_t b = 0; b < alphabet_size; b++)
        numerator[k] += f_powers[b] * log(f[b]);
      //--------
      numerator[k] *= (gamma / (sequences[i].length() - M.size() + 1));
      //--------
    }
    //--------
    numerator.push_back(has_motif * (1 - gamma));
    //--------
    double denominator = smithlab::log_sum_log_vec(numerator, numerator.size());
    for (size_t k = 0; k < indicators[i].size(); k++)
      indicators[i][k] = exp(numerator[k] - denominator);
  }
}

void Model::maximization_seq(const vector<string> &sequences,
    const vector<vector<double> > &indicators, vector<double> &zoops_i) {

  const size_t n = sequences.size();
  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0.0));

  size_t total_bg_length = 0;
  //--------
  double new_gamma = 0.0;
  //--------
  for (size_t i = 0; i < n; i++) {
    size_t l = sequences[i].length();
    total_bg_length += (indicators[i].size() - 1);
    //--------
    double new_Q = 0.0;
    //--------
    for (size_t k = 0; k < indicators[i].size(); k++) {
      //--------
      new_Q += indicators[i][k];
      //--------
      for (size_t j = 0; j < l; j++)
        if (j >= k && j < k + M.size())
          nb[j - k + 1][base2int(sequences[i][j])] += indicators[i][k];
        else
          nb[0][base2int(sequences[i][j])] += indicators[i][k];
    }
    //--------
    zoops_i[i] = new_Q;
    new_gamma += zoops_i[i];
    //--------
  }
  //--------
  gamma = max(new_gamma / n, std::numeric_limits<double>::min());
  //--------

  for (size_t b = 0; b < alphabet_size; b++) {
    f[b] = max(nb[0][b] / total_bg_length, std::numeric_limits<double>::min());
    for (size_t j = 0; j < M.size(); j++)
      M[j][b] = max(nb[j + 1][b] / n, std::numeric_limits<double>::min());
  }
}

void Model::set_model(const string &motif) {

  f.resize(alphabet_size, 1.0 / alphabet_size);
  for (size_t i = 0; i < motif.length(); ++i) {
    M[i].resize(alphabet_size, 0.17);
    M[i][base2int_RNA(motif[i])] = 0.49;
  }
}

void Model::determineStartingPoint_best_kmer(const vector<string> &sequences,
    vector<kmer_info> &top_five_kmers) {

  size_t n_top_kmers = 5;
  size_t k_value = M.size();

  const size_t n_kmers = (1ul << 2 * k_value);
  cerr << "N SEQUENCES:\t" << sequences.size() << endl
 << "K MERS TO CHECK:\t" << n_kmers << endl;

  vector<double> base_comp;
  compute_base_comp(sequences, base_comp);

  cerr << "BASE COMP" << endl;
  for (size_t i = 0; i < base_comp.size(); ++i)
cerr << int2base(i) << '\t' << base_comp[i] << endl;

  vector<size_t> lengths;
  for (size_t i = 0; i < sequences.size(); ++i)
    lengths.push_back(sequences[i].length());

  std::priority_queue<kmer_info, vector<kmer_info>, std::greater<kmer_info> > best_kmers;

  cerr << "EVALUATING K-MERS" << endl;

  for (size_t i = 0; i < n_kmers; ++i) {
    const string kmer(i2mer(k_value, i));
    const double expected = expected_seqs_with_kmer(kmer, base_comp, lengths);
    const size_t observed = count_seqs_with_kmer(kmer, sequences);
    best_kmers.push(kmer_info(kmer_info(kmer, expected, observed)));
    if (best_kmers.size() > n_top_kmers)
      best_kmers.pop();
  }

  while (!best_kmers.empty()) {
    top_five_kmers.push_back(best_kmers.top());
    best_kmers.pop();
  }
  reverse(top_five_kmers.begin(), top_five_kmers.end());
}

double Model::compute_kmer_prob(const string &kmer,
    const vector<double> &base_comp) {
  double prob = 1.0;
  for (size_t i = 0; i < kmer.length(); ++i)
    prob *= base_comp[base2int(kmer[i])];
  return prob;
}

// this is just Poisson probability for 0 observations
double Model::prob_no_occurrence(const double prob, const size_t seq_len) {
  return std::exp(-static_cast<int>(seq_len) * prob);
}

double Model::expected_seqs_with_kmer(const string &kmer,
    const vector<double> &base_comp, const vector<size_t> &lengths) {
  const double p = compute_kmer_prob(kmer, base_comp);
  double expected = 0.0;
  for (size_t i = 0; i < lengths.size(); ++i)
    expected += (1.0 - prob_no_occurrence(p, lengths[i]));
  return expected;
}

size_t Model::count_seqs_with_kmer(const string &kmer,
    const vector<string> &sequences) {
  size_t count = 0;
  for (size_t i = 0; i < sequences.size(); ++i) {
    bool has_kmer = false;
    const size_t lim = sequences[i].length() - kmer.length() + 1;
    for (size_t j = 0; j < lim && !has_kmer; ++j)
      has_kmer = !sequences[i].compare(j, kmer.length(), kmer);
    count += has_kmer;
  }
  return count;
}

void Model::compute_base_comp(const vector<string> &sequences,
    vector<double> &base_comp) {
  base_comp.resize(smithlab::alphabet_size, 0.0);
  size_t total = 0;
  for (size_t i = 0; i < sequences.size(); ++i) {
    for (size_t j = 0; j < sequences[i].length(); ++j)
      ++base_comp[base2int(sequences[i][j])];
    total += sequences[i].length();
  }
  std::transform(
      base_comp.begin(), base_comp.end(), base_comp.begin(),
      std::bind2nd(std::divides<double>(), total));
}
