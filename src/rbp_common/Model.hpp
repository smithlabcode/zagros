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
  Model(const size_t motif_width, const std::vector<GenomicRegion> &regions);

  void
  init_model(const size_t motif_width, const int delta_param);

  /*** parameter fitting ***/

  double calculateLogL(const std::vector<std::string> &S,
      const std::vector<std::vector<double> > &I);

  double calculateLogL(const std::vector<std::string> &S,
      const std::vector<std::vector<size_t> > &D,
      const std::vector<std::vector<double> > &I);

  double calculateLogL_seq(const std::vector<std::string> &S,
      const std::vector<std::vector<size_t> > &D,
      const std::vector<std::vector<double> > &I);

  void expectation_seq(const std::vector<std::string> &S,
      std::vector<std::vector<double> > &I);

  void maximization_seq(const std::vector<std::string> &S,
      const std::vector<std::vector<double> > &I);

  void expectation_seq_de(const std::vector<std::string> &S,
      const std::vector<GenomicRegion> &regions,
      const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I);

  void maximization_seq_de(const std::vector<std::string> &S,
      const std::vector<GenomicRegion> &regions,
      const std::vector<std::vector<size_t> > &D,
      const std::vector<std::vector<double> > &I);

  void
  expectation_str(const std::vector<std::vector<double> > &T,
      const std::vector<double> &wnSecStr,
      std::vector<std::vector<double> > &I);

  void
  maximization_str(const std::vector<std::vector<double> > &T,
      const std::vector<double> &wnSecStr,
      const std::vector<std::vector<double> > &I);
  void
  expectation_seq_str(const std::vector<std::string> &S,
      const std::vector<std::vector<double> > &T,
      const std::vector<double> &wnSecStr,
      std::vector<std::vector<double> > &I);

  void
  maximization_seq_str(const std::vector<std::string> &S,
      const std::vector<std::vector<double> > &T,
      const std::vector<double> &wnSecStr,
      const std::vector<std::vector<double> > &I);

  void expectation_de(const std::vector<GenomicRegion> &regions,
      const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I);

  void maximization_de(const std::vector<GenomicRegion> &regions,
      const std::vector<std::vector<size_t> > &D,
      const std::vector<std::vector<double> > &I);

  void expectation_maximization(const size_t max_iterations,
      const double tolerance, const std::vector<GenomicRegion> &regions,
      std::vector<std::string> &S, const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I, bool s, bool t, size_t d,
      std::string &file_name_base);

  void expectation_maximization_seq(const size_t max_iterations,
      const double tolerance, std::vector<std::string> &S,
      const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I);

  void expectation_maximization_seq_str(const size_t max_iterations,
      const double tolerance, std::vector<std::string> &S,
      const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I, std::string &file_name_base);

  void expectation_maximization_seq_de(const size_t max_iterations,
      const double tolerance, const std::vector<GenomicRegion> &regions,
      std::vector<std::string> &S, const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I);

  void expectation_maximization_str_de(const size_t max_iterations,
      const double tolerance, const std::vector<GenomicRegion> &regions,
      std::vector<std::string> &S, const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I, std::string &file_name_base);

  void expectation_maximization_seq_str_de(const size_t max_iterations,
      const double tolerance, const std::vector<GenomicRegion> &regions,
      std::vector<std::string> &S, const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I, std::string &file_name_base);

  void find_delta(const std::vector<std::string> &sequences,
      const std::vector<GenomicRegion> &regions,
      const std::vector<std::vector<size_t> > &D);

  void expectation_maximization_seq(const size_t max_iterations,
      const double tolerance, std::vector<std::string> &S,
      std::vector<std::vector<double> > &I);

  void
  sampling(const std::vector<std::string> &S,
      std::vector<std::vector<double> > &I);
  void
  update(const std::vector<std::string> &S,
      const std::vector<std::vector<double> > &I);
  void
  gibbs_sampling(const size_t max_iterations, const double tolerance,
      const std::vector<std::string> &S,
      const std::vector<std::vector<size_t> > &D,
      std::vector<std::vector<double> > &I);

  std::string
  print_model(const std::string &motif_name,
      const std::vector<GenomicRegion> &regions,
      const std::vector<std::string> &sequences,
      const std::vector<std::vector<double> > &indicators);

  void
  prepare_output(std::vector<std::string> &seqs,
      const std::vector<std::vector<double> > &indicators,
      const std::vector<std::vector<size_t> > &diagnostic_events,
      const std::string &base_file);

  std::string make_pwm();

  std::string make_structure_profile();

  std::string make_location_profile(
      const std::vector<std::vector<double> > &indicators);

  std::string make_seqs_profile(
      const std::vector<std::vector<double> > &indicators,
      const std::vector<std::string> &sequences);

  std::string make_de_profile(
      const std::vector<std::vector<double> > &indicators,
      const std::vector<std::vector<size_t> > &diagnostic_events,
      const std::vector<std::string> &sequences);

  void generate_profile(const std::vector<std::vector<double> > &indicators,
      const std::vector<std::vector<size_t> > &diagnostic_events,
      const std::vector<std::string> &seqs, const std::string &base_file,
      const std::string option);

  int get_delta() {
    return delta;
  }

  void set_delta(int param) {
    delta = param;
  }

private:
  /*** Private member constants ***/
  std::vector<std::vector<double> > M;
  std::vector<double> lambda;
  std::vector<double> f;
  double p;
  int delta;

  std::vector<double> structure_profile;
};

#endif
