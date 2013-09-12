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
#include <tr1/unordered_map>

#include "smithlab_utils.hpp"
#include "Model.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;

using std::tr1::unordered_map;

using smithlab::alphabet_size;


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
Model::set_model_uniform(const size_t width, Model &model) {
  model.matrix.clear();
  model.matrix.resize(width, vector<double>(alphabet_size, 1.0/alphabet_size));
  model.lambda = vector<double>(width, 0.5);
  model.f = vector<double>(alphabet_size, 1.0/alphabet_size);
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
	       // This value needs to be changed
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
			    const vector<vector<double> > &site_indic) const{
  
  vector<vector<double> > nb_fg;
  vector<double> nb_bg;
  calculate_number_of_bases_fg_bg(sequences, site_indic, matrix.size(), 
				  nb_fg, nb_bg);
  
  double ret = 0.0;
  for (size_t i = 0; i < alphabet_size; ++i) {
    ret += nb_bg[i] * log(f[i]);
    for (size_t j = 0; j < matrix.size(); ++j)
      ret += nb_fg[j][i] * log(matrix[j][i]);
  }
  return ret;
}

double
Model::calculate_zoops_log_l(const vector<string> &sequences,
			     const vector<vector<double> > &site_indic,
			     const vector<double> &seq_indic) const {
  

  vector<vector<double> > nb_fg;
  vector<double> nb_bg;
  calculate_number_of_bases_fg_bg(sequences, site_indic, matrix.size(), 
				  nb_fg, nb_bg);
  
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
    ret += (1-seq_indic[i]) * has_no_motif;
    ret += (1-seq_indic[i]) * log(1-gamma);
    ret += seq_indic[i] * log(gamma / (sequences[i].length() - matrix.size() + 1));
  }
  //--------

  return ret;
}


void
Model::set_model_by_word(const double pseudocount, 
			 const string &kmer, Model &model) {
  
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
      const double tot = 
	accumulate(model.matrix[i].begin(), model.matrix[i].end(), 0.0);
      transform(model.matrix[i].begin(), model.matrix[i].end(), 
		model.matrix[i].begin(),
		std::bind2nd(std::divides<double>(), tot));
    }
}


void
Model::expectation_maximization(const vector<string> &sequences,
				const vector<vector<size_t> > &diagnostic_events,
				const vector<double> &secondary_structure,
				vector<vector<double> > &site_indic, 
				vector<double> &seq_indic) {
  expectation_maximization_seq(sequences, site_indic, seq_indic);
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
    else f_powers[base]++;
    assert(std::isfinite(f_powers[base]) && std::isfinite(num));
  }
  for (size_t b = 0; b < alphabet_size; b++)
    num += f_powers[b]*log(freqs[b]);
  num += log(gamma/(seq.length() - matrix.size() + 1.0));
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

  const double denominator = 
    smithlab::log_sum_log_vec(numerator, numerator.size());
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
  
  vector<vector<double> > nb_fg(matrix.size(),
				// This value needs to be changed
				vector<double>(alphabet_size, pseudocount));
  vector<double> nb_bg(alphabet_size, pseudocount);
  calculate_number_of_bases_fg_bg(sequences, site_indic, matrix.size(), 
				  nb_fg, nb_bg);
  
  for (size_t i = 0; i < matrix.size(); ++i) {
    const double total = accumulate(nb_fg[i].begin(), nb_fg[i].end(), 0.0);
    transform(nb_fg[i].begin(), nb_fg[i].end(), matrix[i].begin(),
	      std::bind2nd(std::divides<double>(), total));
  }
  
  
  const double total = accumulate(nb_bg.begin(), nb_bg.end(), 0.0);
  transform(nb_bg.begin(), nb_bg.end(), freq.begin(),
	    std::bind2nd(std::divides<double>(), total));
  
  gamma = accumulate(seq_indic.begin(), seq_indic.end(), 0.0)/sequences.size();
}



void 
Model::expectation_maximization_seq(const vector<string> &sequences,
				    vector<vector<double> > &site_indic, 
				    vector<double> &seq_indic) {
  
  double prev_score = std::numeric_limits<double>::max();
  for (size_t i = 0; i < max_iterations; ++i) {
    
    expectation_seq(sequences, matrix, f, gamma, site_indic, seq_indic);
    maximization_seq(sequences, site_indic, seq_indic, matrix, f, gamma);

    const double score = 
      calculate_zoops_log_l(sequences, site_indic, seq_indic);
    
    if ((prev_score - score) / prev_score < tolerance) {
      break;
    }
  }
}
