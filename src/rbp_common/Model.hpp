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

#include <vector>
#include <string>
#include "smithlab_utils.hpp"

struct Model {

  Model() :
      f(
          std::vector<double>(
              smithlab::alphabet_size, 1.0 / smithlab::alphabet_size)) {
  }

  void
  expectation_maximization(const std::vector<std::string> &sequences,
                           const std::vector<std::vector<size_t> > &diagnostic_events,
                           const std::vector<std::vector<double> > &secondary_structure,
                           std::vector<std::vector<double> > &site_indic,
                           std::vector<double> &seq_indic);

  void
  expectation_maximization_seq(const std::vector<std::string> &sequences,
                               std::vector<std::vector<double> > &site_indic,
                               std::vector<double> &seq_indic);

  void
  expectation_maximization_seq_str(const std::vector<std::string> &sequences,
                                   const std::vector<std::vector<double> > &secondary_structure,
                                   std::vector<std::vector<double> > &site_indic,
                                   std::vector<double> &seq_indic);

  double
  calculate_oops_log_l(const std::vector<std::string> &sequences,
                       const std::vector<std::vector<double> > &site_indic) const;

  double
  calculate_zoops_log_l(const std::vector<std::string> &sequences,
                        const std::vector<std::vector<double> > &site_indic,
                        const std::vector<double> &seq_indic) const;

  double
  calculate_zoops_log_l(const std::vector<std::string> &sequences,
                        const std::vector<std::vector<double> > &secondary_structure,
                        const std::vector<std::vector<double> > &site_indic,
                        const std::vector<double> &seq_indic) const;

  size_t size() const {
    return matrix.size();
  }

  /* THESE COULD EASILY BE CONSTRUCTORS, BUT DANGEROUS TO GIVE
   SEMANTICS TO ANY DEFAULT IN THIS CASE */
  static void
  set_model_uniform(const size_t width,
                    Model &model);
  static void
  set_model_by_word(const double pseudo,
                    const std::string &kmer,
                    Model &model);

  // instance variables
  std::vector<std::vector<double> > matrix;
  std::vector<double> motif_sec_str;
  std::vector<double> f;
  double f_sec_str;
  double p;
  double gamma;
  int delta;

  // class constants (probably some should be adjustable)
  static const double pseudocount = 0.1;
  static const double zoops_threshold = 0.5;
  static const size_t max_iterations = 5;
  static const double tolerance = 1e-10;
};

#endif
