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
    useStructure = false;
    useDEs = false;
    opt_delta = true;
    opt_geo = true;
  }

  void
  expectationMax(const std::vector<std::string> &sequences,
                 const std::vector<std::vector<double> > &diagnostic_events,
                 const std::vector<std::vector<std::vector<double> > > &diag_values,
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

  void
  expectation_maximization_seq_de(const std::vector<std::string> &sequences,
                                  const std::vector<std::vector<double> > &diagnostic_events,
                                  const std::vector<std::vector<std::vector<double> > > &diag_values,
                                  std::vector<std::vector<double> > &site_indic,
                                  std::vector<double> &seq_indic,
                                  const bool holdDelta);

  void
  expectationMax_SeqStrDE(const std::vector<std::string> &sequences,
                          const std::vector<std::vector<double> > &secStructure,
                          const std::vector<std::vector<double> > &diagnostic_events,
                          const std::vector<std::vector<std::vector<double> > > &diag_values,
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
                        const std::vector<std::vector<double> > &diagnostic_events,
                        const std::vector<std::vector<std::vector<double> > > &diag_values,
                        const std::vector<std::vector<double> > &site_indic,
                        const std::vector<double> &seq_indic) const;

  double
  calculate_zoops_log_l(const std::vector<std::string> &sequences,
                        const std::vector<std::vector<double> > &secondary_structure,
                        const std::vector<std::vector<double> > &diagnostic_events,
                        const std::vector<std::vector<std::vector<double> > > &diag_values,
                        const std::vector<std::vector<double> > &site_indic,
                        const std::vector<double> &seq_indic) const;

  void
  estimateDelta (const std::vector<std::string> &seqs,
                 const std::vector<std::vector<double> > &diagnostic_events,
                 const std::vector<std::vector<std::vector<double> > > &diag_values);

  std::string
  toString_pwm() const;

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

  // Flags for identifying the extra information, if they are used
  bool useStructure;
  bool useDEs;

  // flags for what to optimise
  bool opt_delta;
  bool opt_geo;

  // class constants (probably some should be adjustable)
  static const double DE_WEIGHT;

  static const size_t max_iterations = 10;
  static const double pseudocount;
  static const double zoops_threshold;
  static const double tolerance;

  static const double DEFAULT_GEO_P;
  static const double MIN_GEO_P;
  static const double MAX_GEO_P;

  static const int MIN_DELTA = -8;
  static const int MAX_DELTA = 8;
  static const int DEFAULT_DELTA = 0;
  static const bool HOLD_DELTA_FIXED = true;

  static const size_t DEBUG_LEVEL = 1;
};

#endif
