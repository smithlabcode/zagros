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
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Model.hpp"
#include "IO.hpp"
#include "RNA_Utils.hpp"

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
}

void Model::init_model(const size_t motif_width, const int delta_param) {

  M.clear();
  for (size_t i = 0; i < motif_width; ++i)
    M.push_back(vector<double>(alphabet_size, 1.0 / alphabet_size));
  lambda.resize(motif_width, 0.5);
  f.resize(alphabet_size, 1.0 / alphabet_size);
  p = 0.5;
  set_delta(delta_param);
}

double Model::calculateLogL(const vector<string> &S,
    const vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0.0));

  for (size_t i = 0; i < n; i++) {
    const size_t l = S[i].length();
    for (size_t k = 0; k < I[i].size(); k++)
      for (size_t j = 0; j < l; j++)
        if (j >= k && j < k + M.size())
          nb[j - k + 1][base2int(S[i][j])] += I[i][k];
        else
          nb[0][base2int(S[i][j])] += I[i][k];
  }

  double ret = 0.0;
  for (size_t b = 0; b < alphabet_size; b++) {
    ret += nb[0][b] * log(f[b]);
    for (size_t j = 0; j < M.size(); j++)
      ret += nb[j + 1][b] * log(M[j][b]);
  }

  return ret;
}

double Model::calculateLogL(const vector<string> &S,
    const vector<vector<size_t> > &D, const vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0.0));

  for (size_t i = 0; i < n; i++) {
    const size_t l = S[i].length();
    for (size_t k = 0; k < I[i].size(); k++)
      for (size_t j = 0; j < l; j++)
        if (j >= k && j < k + M.size())
          nb[j - k + 1][base2int(S[i][j])] += I[i][k];
        else
          nb[0][base2int(S[i][j])] += I[i][k];
  }

  double ret = 0.0;
  for (size_t b = 0; b < alphabet_size; b++) {
    ret += nb[0][b] * log(f[b]);
    for (size_t j = 0; j < M.size(); j++)
      ret += nb[j + 1][b] * log(M[j][b]);
  }

  for (size_t i = 0; i < n; i++) {
    if (D[i].size() > 0) {
      for (size_t k = 0; k < I[i].size(); k++) {
        double power = 0.0;
        for (size_t j = 0; j < D[i].size(); j++)
          power += abs(D[i][j] - (k + delta));
        assert(std::isfinite(power));
        ret += I[i][k] * ((power * log(1 - p)) + (D[i].size() * log(p)));
      }
    }
  }

  return ret;
}

double Model::calculateLogL_seq(const vector<string> &S,
    const vector<vector<size_t> > &D, const vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0.0));

  for (size_t i = 0; i < n; i++) {
    const size_t l = S[i].length();
    for (size_t k = 0; k < I[i].size(); k++)
      for (size_t j = 0; j < l; j++)
        if (j >= k && j < k + M.size())
          nb[j - k + 1][base2int(S[i][j])] += I[i][k];
        else
          nb[0][base2int(S[i][j])] += I[i][k];
  }

  double ret = 0.0;
  for (size_t b = 0; b < alphabet_size; b++) {
    ret += nb[0][b] * log(f[b]);
    for (size_t j = 0; j < M.size(); j++)
      ret += nb[j + 1][b] * log(M[j][b]);
  }

  return ret;
}

void Model::maximization_seq_str(const vector<string> &S,
    const vector<vector<double> > &T, const vector<double> &wnSecStr,
    const vector<vector<double> > &I) {

  const size_t n = S.size();

  double lambda_tmp;
  for (size_t j = 0; j < M.size(); j++) {
    lambda_tmp = 0;
    for (size_t i = 0; i < n; i++) {
      for (size_t k = 0; k < I[i].size(); k++)
        lambda_tmp += I[i][k] * T[i][k + j];
    }
    lambda[j] = lambda_tmp / n;
  }

  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0.0));

  size_t total_bg_length = 0;
  for (size_t i = 0; i < n; i++) {
    size_t l = S[i].length();
    total_bg_length += (I[i].size() - 1);

    for (size_t k = 0; k < I[i].size(); k++)
      for (size_t j = 0; j < l; j++)
        if (j >= k && j < k + M.size())
          nb[j - k + 1][base2int(S[i][j])] += I[i][k];
        else
          nb[0][base2int(S[i][j])] += I[i][k];
  }

  for (size_t b = 0; b < alphabet_size; b++) {
    f[b] = max(nb[0][b] / total_bg_length, std::numeric_limits<double>::min());
    for (size_t j = 0; j < M.size(); j++)
      M[j][b] = max(nb[j + 1][b] / n, std::numeric_limits<double>::min());
  }

}

void Model::expectation_seq_str(const vector<string> &S,
    const vector<vector<double> > &T, const vector<double> &wnSecStr,
    vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<vector<double> > alpha = RNAUtils::calculateUnpairedState(
      S, T, wnSecStr);
  vector<vector<double> > gamma = RNAUtils::calculatePairedState(
      S, T, wnSecStr);

  vector<double> denominator;
  for (size_t i = 0; i < n; i++) {

    const size_t l = S[i].length();
    vector<double> numerator(I[i].size(), 0.0);

    for (size_t k = 0; k < I[i].size(); k++) {

      vector<double> f_powers_ss(alphabet_size, 0.0);
      vector<double> f_powers_ds(alphabet_size, 0.0);
      for (size_t j = 0; j < l; j++) {

        const char base = base2int(S[i][j]);
        if (j >= k && j < k + M.size()) {
          numerator[k] += (alpha[i][j] * log(M[j - k][base] * lambda[j - k]));
          numerator[k] += (gamma[i][j] * log(M[j - k][base] * lambda[j - k]));
        } else {
          f_powers_ss[base2int(S[i][j])] += (alpha[i][j]);
          f_powers_ds[base2int(S[i][j])] += (gamma[i][j]);
        }

        assert(std::isfinite(f_powers_ss[base]));
        if (!std::isfinite(numerator[k])) {
          cerr << "Here - Seq_Str!" << endl;
          exit(1);
        }
      }
      for (size_t b = 0; b < alphabet_size; b++) {
        numerator[k] += f_powers_ss[b] * log(f[b] * 0.5);
        numerator[k] += f_powers_ds[b] * log(f[b] * 0.5);
      }
    }

    denominator.push_back(
        smithlab::log_sum_log_vec(numerator, numerator.size()));
    for (size_t k = 0; k < I[i].size(); k++)
        I[i][k] = exp(numerator[k] - denominator[i]);
  }
}

void Model::expectation_seq_str_de(const vector<string> &S,
    const vector<vector<double> > &T, const vector<double> &wnSecStr,
    const vector<vector<size_t> > &D, vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<vector<double> > alpha = RNAUtils::calculateUnpairedState(
      S, T, wnSecStr);
  vector<vector<double> > gamma = RNAUtils::calculatePairedState(
      S, T, wnSecStr);

  vector<double> denominator;
  for (size_t i = 0; i < n; i++) {

    const size_t l = S[i].length();
    vector<double> numerator(I[i].size(), 0.0);

    for (size_t k = 0; k < I[i].size(); k++) {

      vector<double> f_powers_ss(alphabet_size, 0.0);
      vector<double> f_powers_ds(alphabet_size, 0.0);
      for (size_t j = 0; j < l; j++) {

        const char base = base2int(S[i][j]);
        if (j >= k && j < k + M.size()) {
          numerator[k] += (alpha[i][j] * log(M[j - k][base] * lambda[j - k]));
          numerator[k] += (gamma[i][j] * log(M[j - k][base] * lambda[j - k]));
        } else {
          f_powers_ss[base2int(S[i][j])] += (alpha[i][j]);
          f_powers_ds[base2int(S[i][j])] += (gamma[i][j]);
        }

        assert(std::isfinite(f_powers_ss[base]));
        if (!std::isfinite(numerator[k])) {
          cerr << "Here - Seq_Str!" << endl;
          exit(1);
        }
      }
      for (size_t b = 0; b < alphabet_size; b++) {
        numerator[k] += f_powers_ss[b] * log(f[b] * 0.5);
        numerator[k] += f_powers_ds[b] * log(f[b] * 0.5);
      }

      if (D[i].size() > 0) {
        double power = 0.0;
        for (size_t j = 0; j < D[i].size(); j++)
          power += abs(D[i][j] - (k + delta));
        assert(std::isfinite(power));
        if (!std::isfinite(numerator[k])) {
          cerr << "DE __  Here!" << endl;
          exit(1);
        }
        numerator[k] = (numerator[k]
            + ((power * log(1 - p)) + (D[i].size() * log(p))));
      }
    }

    denominator.push_back(
        smithlab::log_sum_log_vec(numerator, numerator.size()));
    for (size_t k = 0; k < I[i].size(); k++)
        I[i][k] = exp(numerator[k] - denominator[i]);
  }
}

void Model::expectation_str_de(const vector<vector<double> > &T,
    const vector<double> &wnSecStr,
    const vector<vector<size_t> > &D,
    vector<vector<double> > &I) {

}

void Model::maximization_str(const vector<vector<double> > &T,
    const vector<double> &wnSecStr, const vector<vector<double> > &I) {

}

void Model::expectation_seq(const vector<string> &S,
    vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<double> denominator;

  for (size_t i = 0; i < n; i++) {

    const size_t l = S[i].length();
    vector<double> numerator(I[i].size(), 0.0);

    for (size_t k = 0; k < I[i].size(); k++) {

      vector<double> f_powers(alphabet_size, 0.0);
      for (size_t j = 0; j < l; j++) {

        const char base = base2int(S[i][j]);
        if (j >= k && j < k + M.size()) {
          numerator[k] += log(M[j - k][base]);
        } else
          f_powers[base2int(S[i][j])]++;

        assert(std::isfinite(f_powers[base]));
        if (!std::isfinite(numerator[k])) {
          cerr << "Here!" << endl;
          exit(1);
        }
      }
      for (size_t b = 0; b < alphabet_size; b++)
        numerator[k] += f_powers[b] * log(f[b]);
    }

    denominator.push_back(
        smithlab::log_sum_log_vec(numerator, numerator.size()));
    for (size_t k = 0; k < I[i].size(); k++)
      I[i][k] = exp(numerator[k] - denominator[i]);
  }
}

void Model::maximization_seq(const vector<string> &S,
    const vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0.0));

  size_t total_bg_length = 0;
  for (size_t i = 0; i < n; i++) {
    size_t l = S[i].length();
    total_bg_length += (I[i].size() - 1);

    for (size_t k = 0; k < I[i].size(); k++)
      for (size_t j = 0; j < l; j++)
        if (j >= k && j < k + M.size())
          nb[j - k + 1][base2int(S[i][j])] += I[i][k];
        else
          nb[0][base2int(S[i][j])] += I[i][k];
  }

  for (size_t b = 0; b < alphabet_size; b++) {
    f[b] = max(nb[0][b] / total_bg_length, std::numeric_limits<double>::min());
    for (size_t j = 0; j < M.size(); j++)
      M[j][b] = max(nb[j + 1][b] / n, std::numeric_limits<double>::min());
  }
}

void Model::maximization_de(const vector<GenomicRegion> &regions,
    const vector<vector<size_t> > &D, const vector<vector<double> > &I) {
  const size_t n = I.size();
  double total_sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    if (D[i].size() > 0) {
      double indicators_sum = 0.0;
      for (size_t k = 0; k < I[i].size(); k++) {
        double de_sum = 0.0;
        for (size_t j = 0; j < D[i].size(); j++)
          de_sum += abs(D[i][j] - (k + delta));
        indicators_sum += (I[i][k] * D[i].size()) / (D[i].size() + de_sum);
      }
      total_sum += indicators_sum;
    }
  }
  p = max(total_sum / n, std::numeric_limits<double>::min());
}

void Model::expectation_seq_de(const vector<string> &S,
    const vector<GenomicRegion> &regions, const vector<vector<size_t> > &D,
    vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<double> denominator;
  for (size_t i = 0; i < n; i++) {

    const size_t l = S[i].length();
    vector<double> numerator(I[i].size(), 0.0);

    for (size_t k = 0; k < I[i].size(); k++) {
      vector<double> f_powers(alphabet_size, 0.0);
      for (size_t j = 0; j < l; j++) {
        const char base = base2int(S[i][j]);
        if (j >= k && j < k + M.size()) {
          numerator[k] += log(M[j - k][base]);
        } else
          f_powers[base2int(S[i][j])]++;
        assert(std::isfinite(f_powers[base]));
        if (!std::isfinite(numerator[k])) {
          cerr << "Here!" << endl;
          exit(1);
        }
      }
      for (size_t b = 0; b < alphabet_size; b++)
        numerator[k] += f_powers[b] * log(f[b]);

      if (D[i].size() > 0) {
        double power = 0.0;
        for (size_t j = 0; j < D[i].size(); j++)
          power += abs(D[i][j] - (k + delta));
        assert(std::isfinite(power));
        if (!std::isfinite(numerator[k])) {
          cerr << "DE __  Here!" << endl;
          exit(1);
        }
        numerator[k] = (numerator[k]
            + ((power * log(1 - p)) + (D[i].size() * log(p))));
      }
    }
    denominator.push_back(
        smithlab::log_sum_log_vec(numerator, numerator.size()));
    for (size_t k = 0; k < I[i].size(); k++)
      I[i][k] = exp(numerator[k] - denominator[i]);
  }
}

void Model::expectation_maximization(const size_t max_iterations,
    const double tolerance, const vector<GenomicRegion> &regions,
    vector<string> &S, const vector<vector<size_t> > &D,
    vector<vector<double> > &I, bool s, bool t, size_t d,
    string &file_name_base) {

  cerr << "Fitting the shifting parameter..." << endl;
  find_delta(S, regions, D);

  if (s && t && !d)
    expectation_maximization_seq_str(
        max_iterations, tolerance, S, D, I, file_name_base);
  else if (s && !t && d) {
    expectation_maximization_seq_de(
        max_iterations, tolerance, regions, S, D, I);
  } else if (!s && t && d) {
    expectation_maximization_str_de(
        max_iterations, tolerance, regions, S, D, I, file_name_base);
  } else if (s && t && d) {
    expectation_maximization_seq_str_de(
        max_iterations, tolerance, regions, S, D, I, file_name_base);
  } else
    expectation_maximization_seq(max_iterations, tolerance, S, D, I);
}

void Model::find_delta(const vector<string> &sequences,
    const vector<GenomicRegion> &regions, const vector<vector<size_t> > &D) {

  vector<double> ll_delta;
  for (int delta_param = -10; delta_param < 11; ++delta_param) {
    cerr << "\r"
        << int(
            100 * ll_delta.size()
                / 21)
        << "% completed..." << std::flush;
    init_model(M.size(), delta_param);
    const size_t max_iterations = 10;
    const double tolerance = 1e-10;
    vector<vector<double> > indicators;
    for (size_t i = 0; i < sequences.size(); ++i) {
      const size_t n_pos = sequences[i].length() - M.size() + 1;
      indicators.push_back(vector<double>(n_pos, 1 / n_pos));
    }
    double prev_score = std::numeric_limits<double>::max();
    for (size_t i = 0; i < max_iterations; ++i) {
      expectation_seq_de(sequences, regions, D, indicators);
      maximization_seq(sequences, indicators);
      maximization_de(regions, D, indicators);
      const double score = calculateLogL(sequences, D, indicators);

      if ((prev_score - score) / prev_score < tolerance) {
        break;
      }
      prev_score = score;
    }
    ll_delta.push_back(prev_score);
  }
  cerr << "\r" << "100% completed..." << endl;

  double max_ll = -100000000;
  int max_i = 0;
  for (size_t i = 0; i < ll_delta.size(); i++) {
    if (ll_delta[i] > max_ll) {
      max_ll = ll_delta[i];
      max_i = i - 10;
    }
  }
  init_model(M.size(), max_i);
}

void Model::expectation_maximization_seq(const size_t max_iterations,
    const double tolerance, vector<string> &S, const vector<vector<size_t> > &D,
    vector<vector<double> > &I) {

  cerr << "Fitting the full model started...";

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq(S, I);
    maximization_seq(S, I);

    const double score = calculateLogL(S, I);

    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
    prev_score = score;
  }

  cerr << "done!" << endl;

  double max_X = -1;
  int max_i = -1;

  for (size_t n = 0; n < S.size(); n++) {
    max_X = -1;
    max_i = -1;
    for (size_t i = 0; i < I[n].size(); i++) {
      if (I[n][i] > max_X) {
        max_X = I[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < M.size(); ++j) {
      M[j][RNAUtils::base2int_RNA(S[n][max_i + j])] += 1;
    }
  }
}

void Model::expectation_maximization_seq_str(const size_t max_iterations,
    const double tolerance, vector<string> &S, const vector<vector<size_t> > &D,
    vector<vector<double> > &I, string &file_name_base) {

//  determineStartingStructure(S);

  vector<vector<double> > fullStrVectorTable;
  vector<vector<vector<double> > > fullStrMatrixTable;
  vector<double> wnSecStr;

  if (IO::str_file_checks_out(S.front(), file_name_base)) {
    IO::fillTables(
        file_name_base, S, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
    IO::trimTables(S, file_name_base, fullStrVectorTable, fullStrMatrixTable);
  } else {
    IO::fillTables(S, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
    IO::trimTables(S, file_name_base, fullStrVectorTable, fullStrMatrixTable);
    IO::saveTables(
        S, file_name_base, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
  }

  cerr << "Fitting the full model started...";

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_str(S, fullStrVectorTable, wnSecStr, I);
    maximization_seq_str(S, fullStrVectorTable, wnSecStr, I);

    const double score = calculateLogL(S, I);

    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
    prev_score = score;
  }

  cerr << "done!" << endl;

  double max_X = -1;
  int max_i = -1;

  for (size_t n = 0; n < S.size(); n++) {
    max_X = -1;
    max_i = -1;
    for (size_t i = 0; i < I[n].size(); i++) {
      if (I[n][i] > max_X) {
        max_X = I[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < M.size(); ++j) {
      M[j][RNAUtils::base2int_RNA(S[n][max_i + j])] += 1;
    }
  }
}

void Model::expectation_maximization_seq_de(const size_t max_iterations,
    const double tolerance, const vector<GenomicRegion> &regions,
    vector<string> &S, const vector<vector<size_t> > &D,
    vector<vector<double> > &I) {

  cerr << "Fitting the full model started...";

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_de(S, regions, D, I);
    maximization_seq(S, I);
    maximization_de(regions, D, I);

    const double score = calculateLogL(S, D, I);

    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
    prev_score = score;
  }
  cerr << "done!" << endl;

  double max_X = -1;
  size_t max_i = 0;

  for (size_t n = 0; n < S.size(); n++) {
    max_X = -1;
    max_i = 0;
    for (size_t i = 0; i < I[n].size(); i++) {
      if (I[n][i] > max_X) {
        max_X = I[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < M.size(); ++j) {
      M[j][RNAUtils::base2int_RNA(S[n][max_i + j])] += 1;
    }
  }
}

void Model::expectation_maximization_str_de(const size_t max_iterations,
    const double tolerance, const vector<GenomicRegion> &regions,
    vector<string> &S, const vector<vector<size_t> > &D,
    vector<vector<double> > &I, string &file_name_base) {

  vector<vector<double> > fullStrVectorTable;
  vector<vector<vector<double> > > fullStrMatrixTable;
  vector<double> wnSecStr;

  if (IO::str_file_checks_out(S.front(), file_name_base)) {
    IO::fillTables(
        file_name_base, S, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
    IO::trimTables(S, file_name_base, fullStrVectorTable, fullStrMatrixTable);
  } else {
    IO::fillTables(S, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
    IO::trimTables(S, file_name_base, fullStrVectorTable, fullStrMatrixTable);
    IO::saveTables(
        S, file_name_base, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
  }

  cerr << "Fitting the full model started...";

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_str_de(fullStrVectorTable, wnSecStr, D, I);
    maximization_str(fullStrVectorTable, wnSecStr, I);
    maximization_de(regions, D, I);

    const double score = calculateLogL(S, D, I);

    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
    prev_score = score;
  }

  cerr << "done!" << endl;

  double max_X = -1;
  int max_i = -1;

  for (size_t n = 0; n < S.size(); n++) {
    max_X = -1;
    max_i = -1;
    for (size_t i = 0; i < I[n].size(); i++) {
      if (I[n][i] > max_X) {
        max_X = I[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < M.size(); ++j) {
      M[j][RNAUtils::base2int_RNA(S[n][max_i + j])] += 1;
    }
  }
}

void Model::expectation_maximization_seq_str_de(const size_t max_iterations,
    const double tolerance, const vector<GenomicRegion> &regions,
    vector<string> &S, const vector<vector<size_t> > &D,
    vector<vector<double> > &I, string &file_name_base) {

  vector<vector<double> > fullStrVectorTable;
  vector<vector<vector<double> > > fullStrMatrixTable;
  vector<double> wnSecStr;

  if (IO::str_file_checks_out(S.front(), file_name_base)) {
    IO::fillTables(
        file_name_base, S, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
    IO::trimTables(S, file_name_base, fullStrVectorTable, fullStrMatrixTable);
  } else {
    IO::fillTables(S, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
    IO::trimTables(S, file_name_base, fullStrVectorTable, fullStrMatrixTable);
    IO::saveTables(
        S, file_name_base, fullStrVectorTable, fullStrMatrixTable, wnSecStr);
  }

  cerr << "Fitting the full model started...";

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    expectation_seq_str_de(S, fullStrVectorTable, wnSecStr, D, I);
    maximization_seq_str(S, fullStrVectorTable, wnSecStr, I);
    maximization_de(regions, D, I);

    const double score = calculateLogL(S, D, I);

    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
    prev_score = score;
  }

  cerr << "done!" << endl;

  double max_X = -1;
  int max_i = -1;

  for (size_t n = 0; n < S.size(); n++) {
    max_X = -1;
    max_i = -1;
    for (size_t i = 0; i < I[n].size(); i++) {
      if (I[n][i] > max_X) {
        max_X = I[n][i];
        max_i = i;
      }
    }
    for (size_t j = 0; j < M.size(); ++j) {
      M[j][RNAUtils::base2int_RNA(S[n][max_i + j])] += 1;
    }
  }
}

void Model::sampling(const vector<string> &S, vector<vector<double> > &I) {

  const size_t n = S.size();
  vector<double> denominator(n, 0.0);

  vector<size_t> indices;
  for (size_t i = 0; i < n; i++)
    indices.push_back(i);

  srand(time(0) + getpid());

  for (size_t d = 0; d < n; d++) {

    const size_t r = rand() % indices.size();
    const size_t idx_r = indices[r];
    vector<double> numerator(I[idx_r].size(), 0.0);
    double m = std::numeric_limits<double>::min();

    vector<string> S_nc = S;
    S_nc.erase(S_nc.begin() + idx_r);
    vector<vector<double> > I_nc = I;
    I_nc.erase(I_nc.begin() + idx_r);

    vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0));
    size_t total_bg_length = 0;
    for (size_t i = 0; i < n - 1; i++) {
      const size_t l = S_nc[i].length();
      total_bg_length += I_nc[i].size() - 1;

      for (size_t k = 0; k < I_nc[i].size(); k++)
        for (size_t j = 0; j < l; j++)
          if (j >= k && j < k + M.size())
            nb[j - k + 1][base2int(S_nc[i][j])] += I_nc[i][k];
          else
            nb[0][base2int(S_nc[i][j])] += I_nc[i][k];
    }

    vector<vector<double> > M_hat(
        M.size(), vector<double>(alphabet_size, 1.0 / alphabet_size));

    vector<double> f_hat(alphabet_size, 1.0 / alphabet_size);
    for (size_t k = 0; k < I[idx_r].size(); k++) {
      for (size_t j = 0; j < M.size(); j++) {
        const char idx = base2int(S[idx_r][k + j]);
        f_hat[idx] = std::max(
            nb[0][idx] / total_bg_length, std::numeric_limits<double>::min());
        M_hat[j][idx] = max(
            nb[j + 1][idx] / (n - 1.0), std::numeric_limits<double>::min());
        numerator[k] += log(M_hat[j][idx] / f_hat[idx]);
      }
      if (numerator[k] > m)
        m = numerator[k];
    }

    for (size_t k = 0; k < I[idx_r].size(); k++)
      if (numerator[k] - m >= -10)
        denominator[idx_r] += exp(numerator[k] - m);

    for (size_t k = 0; k < I[idx_r].size(); k++) {
      if (numerator[k] - m >= -10)
        I[idx_r][k] = exp(numerator[k] - m) / denominator[idx_r];
      else
        I[idx_r][k] = 0.0000005;
    }

    indices.erase(indices.begin() + r);
  }
}

void Model::update(const vector<string> &S, const vector<vector<double> > &I) {

  const size_t n = S.size();

  vector<vector<double> > nb(M.size() + 1, vector<double>(alphabet_size, 0));
  size_t total_bg_length = 0;
  for (size_t i = 0; i < n; i++) {
    const size_t l = S[i].length();
    total_bg_length += I[i].size() - 1;

    double max_X = -1;
    int max_i = -1;
    for (size_t iindx = 0; iindx < I[i].size(); iindx++) {
      if (I[i][iindx] > max_X) {
        max_X = I[i][iindx];
        max_i = iindx;
      }
    }
    size_t k = max_i;
    for (size_t j = 0; j < l; j++)
      if (j >= k && j < k + M.size())
        nb[j - k + 1][base2int(S[i][j])] += 1;
      else
        nb[0][base2int(S[i][j])] += 1;
  }

  for (size_t b = 0; b < alphabet_size; ++b) {
    f[b] = max(nb[0][b] / total_bg_length, std::numeric_limits<double>::min());
    for (size_t j = 0; j < M.size(); j++)
      M[j][b] = max(nb[j + 1][b] / n, std::numeric_limits<double>::min());
  }
}

void Model::gibbs_sampling(const size_t max_iterations, const double tolerance,
    const vector<string> &S, const vector<vector<size_t> > &D,
    vector<vector<double> > &I) {

  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    sampling(S, I);
    const double score = calculateLogL(S, I);
    update(S, I);
    if ((prev_score - score) / prev_score < tolerance)
      break;
    prev_score = score;
  }
  update(S, I);
}

string Model::print_model(const string &motif_name,
    const vector<GenomicRegion> &regions, const vector<string> &sequences,
    const vector<vector<double> > &indicators) {

  stringstream ss;

  const size_t N = sequences.size();
  if (N <= 0) {
    stringstream ss;
    ss << "Building motif alignment failed. Reason: sequences vector is empty";
    throw SMITHLABException(ss.str());
  }

  if (N != indicators.size()) {
    stringstream ss;
    ss << "Building motif alignment failed. Reason: expected " << N
        << "indicator vectors, got only " << indicators.size();
    throw SMITHLABException(ss.str());
  }

  ss << "AC\t" << motif_name << endl;
  ss << "XX" << endl;
  ss << "TY\tMotif" << endl;
  ss << "XX" << endl;
  ss << "P0\tA\tC\tG\tT" << endl;

  for (size_t j = 0; j < M.size(); j++) {
    ss << "0" << j + 1 << "\t";
    for (size_t b = 0; b < alphabet_size - 1; b++)
      ss << (int) (M[j][b]) << "\t";
    ss << (int) (M[j][alphabet_size - 1]) << endl;
  }
  ss << "XX" << endl;
  ss << "AT\tGEO_P=" << p << endl;
  ss << "XX" << endl;

  double max_X = -1;
  size_t max_i = 0;

  for (size_t n = 0; n < N; n++) {
    max_X = -1;
    max_i = 0;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    ss << "BS\t" << sequences[n].substr(max_i, M.size()) << "; "
        << regions[n].get_chrom() << ":" << regions[n].get_start() << "-"
        << regions[n].get_end() << "; " << max_i + 1 << "; " << M.size()
        << ";  ;" << ((regions[n].get_strand() == '+') ? "p; " : "n;") << endl;
  }
  ss << "XX" << endl;
  ss << "//" << endl;

  return ss.str();
}

void Model::prepare_output(vector<string> &seqs,
    const vector<vector<double> > &indicators,
    const vector<vector<size_t> > &diagnostic_events, const string &base_file) {

  string file_name = base_file + "_Rversion.tmp";
  string command = "R --version > " + file_name;
  system(command.c_str());
  file_name = base_file + "_Rversion.tmp";
  ifstream in(file_name.c_str());
  string s;
  getline(in, s);
  in.close();
  command = "rm -f " + file_name;
  system(command.c_str());

  if (s.substr(0, 9) != "R version") {
    cerr << "R is not installed on this system - disabling graphic utilities..."
        << endl;
  } else {
    generate_profile(indicators, diagnostic_events, seqs, base_file, "motif");
    generate_profile(
        indicators, diagnostic_events, seqs, base_file, "locations");
    generate_profile(
        indicators, diagnostic_events, seqs, base_file, "sequences");
    generate_profile(indicators, diagnostic_events, seqs, base_file, "des");
    /*    size_t N = seqs.size();
     vector<vector<double> > fullStrVectorTable;
     IO::fillTables(seqs, fullStrVectorTable);
     IO::trimTables(seqs, fullStrVectorTable);
     for (size_t i = 0; i < seqs.front().length(); ++i)
     structure_profile.push_back(0);
     for (size_t i = 0; i < N; i++) {
     size_t l = seqs[i].length();
     for (size_t j = 0; j < l; ++j)
     structure_profile[j] += fullStrVectorTable[i][j];
     }
     for (size_t i = 0; i < seqs.front().length(); ++i)
     structure_profile[i] = structure_profile[i] / N;
     generate_profile(indicators, diagnostic_events, seqs, base_file, "structure");*/
  }
}

void Model::generate_profile(const vector<vector<double> > &indicators,
    const vector<vector<size_t> > &diagnostic_events,
    const vector<string> &seqs, const string &base_file, const string option) {

  string file_name = base_file + "_" + option + ".dat";
  ofstream outf(file_name.c_str());
  string s;
  if (option == "motif")
    s = make_pwm();
  else if (option == "structure")
    s = make_structure_profile();
  else if (option == "locations")
    s = make_location_profile(indicators);
  else if (option == "sequences")
    s = make_seqs_profile(indicators, seqs);
  else
    s = make_de_profile(indicators, diagnostic_events, seqs);
  outf << s << endl;
  outf.close();
  string command = "cat src/utils/generate_" + option
      + "_profile.R | R --vanilla " + file_name + " ft "
      + convertSizet(seqs.size()) + " " + convertSizet(seqs.front().length())
      + " > " + file_name + "_" + option + ".tmp";
  system(command.c_str());
  command = "rm -f " + file_name;
  system(command.c_str());
  command = "rm -f " + file_name + "_" + option + ".tmp";
  system(command.c_str());
}

string Model::make_pwm() {

  stringstream ss;
  vector<double> sum;
  for (size_t j = 0; j < M.size(); j++) {
    sum.push_back(0);
    for (size_t b = 0; b < alphabet_size; b++)
      sum[j] += M[j][b];
  }

  for (size_t b = 0; b < alphabet_size; b++) {
    for (size_t j = 0; j < M.size() - 1; j++)
      ss << M[j][b] / sum[j] << ",";
    ss << M[M.size() - 1][b] / sum[M.size() - 1] << endl;
  }
  return ss.str();
}

string Model::make_structure_profile() {
  stringstream ss;
  for (size_t i = 0; i < structure_profile.size(); i++)
    ss << structure_profile[i] << endl;
  return ss.str();
}

string Model::make_location_profile(const vector<vector<double> > &indicators) {
  stringstream ss;

  double max_X = -1;
  int max_i = -1;
  for (size_t n = 0; n < indicators.size(); n++) {
    max_X = -1;
    max_i = -1;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    ss << max_i + 1 << endl;
  }
  return ss.str();
}

string Model::make_seqs_profile(const vector<vector<double> > &indicators,
    const vector<string> &sequences) {
  stringstream ss;

  double max_X = -1;
  size_t max_i = 0;
  for (size_t n = 0; n < indicators.size(); n++) {
    max_X = -1;
    max_i = 0;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    for (size_t i = 0; i < sequences[n].length(); i++) {
      string color;
      if (sequences[n].substr(i, 1) == "A")
        color = "green";
      else if (sequences[n].substr(i, 1) == "C")
        color = "blue";
      else if (sequences[n].substr(i, 1) == "G")
        color = "orange";
      else
        color = "red";
      if (i < max_i || (i - max_i) >= M.size())
        ss << 1 << "\t" << sequences[n].substr(i, 1) << "\t" << color << endl;
      else
        ss << 0 << "\t" << sequences[n].substr(i, 1) << "\t" << color << endl;
    }
  }
  return ss.str();
}

string Model::make_de_profile(const vector<vector<double> > &indicators,
    const vector<vector<size_t> > &diagnostic_events,
    const vector<string> &sequences) {
  stringstream ss;

  double max_X = -1;
  size_t max_i = 0;
  for (size_t n = 0; n < indicators.size(); n++) {
    max_X = -1;
    max_i = 0;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }

    vector<double> values(sequences[n].length(), 1.0);

    for (size_t i = 0; i < sequences[n].length(); i++) {
      for (size_t j = 0; j < diagnostic_events[n].size(); ++j)
        if (i == diagnostic_events[n][j])
          values[i]++;
    }

    for (size_t i = 0; i < sequences[n].length(); i++) {
      ss << 1 / values[i] << "\t" << sequences[n].substr(i, 1) << "\t"
          << "yellow" << endl;
    }
  }
  return ss.str();
}

