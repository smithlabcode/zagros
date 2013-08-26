#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <iomanip>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "RNG.hpp"
#include "OptionParser.hpp"

using std::numeric_limits;
using std::string;
using std::vector;
using std::ostream;
using std::stringstream;
using std::istringstream;
using std::endl;
using std::cerr;
using std::cout;

using smithlab::alphabet_size;

static void read_position_weight_matrix(const string &filename,
    vector<vector<double> > &pwm) {

  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("bad PWM file: " + filename);

  string buffer;
  while (getline(in, buffer)) {
    pwm.push_back(vector<double>());
    std::istringstream is(buffer);
    double tmp_val = 0.0;
    while (is >> tmp_val)
      pwm.back().push_back(tmp_val);
  }

  assert(!pwm.empty());
  for (size_t i = 0; i < pwm.size(); ++i) {
    assert(pwm.front().size() == pwm[i].size());
//    const double total = std::accumulate(pwm[i].begin(), pwm[i].end(), 0.0);
//    assert(total == 1.0);
  }
}

static string pwm2dme(const vector<vector<double> > &M, size_t no_seqs) {

  // Check that the parameters make sense
  stringstream ss;

  ss << "AC\tX_MODEL" << endl;
  ss << "XX" << endl;
  ss << "TY\tMotif" << endl;
  ss << "XX" << endl;
  ss << "P0\tA\tC\tG\tT" << endl;

  for (size_t j = 0; j < M.size(); j++) {
    ss << "0" << j + 1 << "\t";
    for (size_t b = 0; b < alphabet_size - 1; b++)
      ss << (int) (M[j][b] * no_seqs) << "\t";
    ss << (int) (M[j][alphabet_size - 1] * no_seqs) << endl;
  }
  ss << "XX" << endl;
  ss << "//" << endl;

  return ss.str();
}

int main(int argc, const char **argv) {

  try {

    size_t motif_width = 0;
    int random_number_seed = numeric_limits<int>::max();
    bool VERBOSE = false;
    string pwm_filename = "NA";

    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse("simrna", "simulate rna for motif discovery"
        "<pwm> <fasta> <bed> <output>");
    opt_parse.add_opt("width", 'w', "motif width", true, motif_width);
    opt_parse.add_opt(
        "load", 'p', "Name of PWM file (default: random generation)", true,
        pwm_filename);
    opt_parse.add_opt(
        "seed", 'S', "random number seed", false, random_number_seed);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (motif_width == 0) {
      cerr << "motif width must be positive" << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }

    const string input_fasta_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> names, sequences;
    read_fasta_file(input_fasta_file, names, sequences);

    vector<vector<double> > pwm;

    read_position_weight_matrix(pwm_filename, pwm);

    vector<vector<double> > motif = pwm;

    string dout = pwm2dme(motif, sequences.size());
    cout << dout;
    // Output the diagnostic events BED file and the updated sequences

  } catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  } catch (SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  } catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
