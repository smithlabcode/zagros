/****
 Description: TODO

 --------------------

 Copyright (C) 2012
 University of Southern California,
 Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

 Authors: Emad Bahrami Samani, Philip J. Uren

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

 --------------------

 Known Bugs:    None

 Revision
 History:       None
 ****/

#ifndef __MODEL_H_
#define __MODEL_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

class Model {
public:
  /*** Constructors, destructors and object initialization ***/
  Model(const size_t motif_width);

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

  void expectation_maximization(const std::vector<std::string> &sequences,
      const std::vector<std::vector<size_t> > &diagnostic_events,
      const std::vector<double> &secondary_structure,
      std::vector<std::vector<double> > &indicators,
      std::vector<double> &zoops_i, const std::string t, const bool d,
      const std::string &file_name_base, const size_t max_iterations,
      const double tolerance);

  int get_delta() {
    return delta;
  }

  void set_delta(int param) {
    delta = param;
  }

  int getDelta() const {
    return delta;
  }

  void setDelta(int delta) {
    this->delta = delta;
  }

  const std::vector<double>& getF() const {
    return f;
  }

  void setF(const std::vector<double>& f) {
    this->f = f;
  }

  double getGamma() const {
    return gamma;
  }

  void setGamma(double gamma) {
    this->gamma = gamma;
  }

  const std::vector<double>& getLambda() const {
    return lambda;
  }

  void setLambda(const std::vector<double>& lambda) {
    this->lambda = lambda;
  }

  const std::vector<std::vector<double> >& getM() const {
    return M;
  }

  double getM_element(const size_t i, const size_t j) const {
    return M[i][j];
  }

  void setM_element(size_t i, size_t j, double val) {
    this->M[i][j] = val;
  }

  size_t get_model_size() const {
    return M.size();
  }

  void setM(const std::vector<std::vector<double> >& m) {
    M = m;
  }

  double getP() const {
    return p;
  }

  void setP(double p) {
    this->p = p;
  }

private:

  /*** Private member constants ***/
  std::vector<std::vector<double> > M;
  std::vector<double> lambda;
  std::vector<double> f;
  double p;
  int delta;

  double gamma;

  struct kmer_info {
    std::string kmer;
    double expected;
    size_t observed;
    kmer_info(const std::string &km, const double ex, const double ob) :
        kmer(km), expected(ex), observed(ob) {
    }
    double score() const {
      return observed / expected;
    }
    bool operator>(const kmer_info &ki) const {
      return score() > ki.score();
    }
  };

  double calculateLogL(const std::vector<std::string> &sequences,
      const std::vector<std::vector<double> > &indicators,
      const std::vector<double> &zoops_i);

  void expectation_maximization_seq(const std::vector<std::string> &sequences,
      std::vector<std::vector<double> > &indicators,
      std::vector<double> &zoops_i, const size_t max_iterations,
      const double tolerance);

  void expectation_seq(const std::vector<std::string> &sequences,
      std::vector<std::vector<double> > &indicators,
      std::vector<double> &zoops_i);

  void maximization_seq(const std::vector<std::string> &sequences,
      const std::vector<std::vector<double> > &indicators,
      std::vector<double> &zoops_i);

  void
  determineStartingPoint_best_kmer(const std::vector<std::string> &sequences,
      std::vector<kmer_info> &top_five_kmers);

  void set_model(const std::string &motif);

  double
  compute_kmer_prob(const std::string &kmer,
      const std::vector<double> &base_comp);

  double
  prob_no_occurrence(const double prob, const size_t seq_len);

  double
  expected_seqs_with_kmer(const std::string &kmer,
      const std::vector<double> &base_comp, const std::vector<size_t> &lengths);

  size_t
  count_seqs_with_kmer(const std::string &kmer,
      const std::vector<std::string> &sequences);

  void
  compute_base_comp(const std::vector<std::string> &sequences,
      std::vector<double> &base_comp);


};

#endif
