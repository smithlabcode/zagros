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

static size_t generate_geometric_random_variable(const Runif &rng,
    const double p) {
  return std::floor(std::log(rng.runif(0.0, 1.0)) / std::log(1.0 - p));
}

static vector<double> generate_dirichlet_random_variables(gsl_rng *gr,
    size_t n) {

  vector<double> gamma;
  vector<double> ret_val;

  double V = 0;
  for (size_t i = 0; i < n; ++i) {
    gamma.push_back(gsl_ran_gamma(gr, 1, 1));
    V += gamma[i];
  }
  for (size_t i = 0; i < n; ++i)
    ret_val.push_back(gamma[i] / V);
  return ret_val;
}

/*static void 
 generate_position_weight_matrix(gsl_rng *gr, size_t n, size_t m, double target_ic, vector<vector<double> > &pwm) {
 double ic = 0;
 double sum = 0;
 vector<double> new_col(m,0.0);

 for (size_t i = 0; pwm.size() < n; ++i) {
 ic = 0;
 new_col = generate_dirichlet_random_variables(gr, m);
 for (size_t j = 0; j < new_col.size(); ++j)
 if (new_col[j] > 0)
 ic += new_col[j] * abs(log(new_col[j])/log(2));
 if (ic > target_ic - 0.001 && ic < target_ic + 0.001) {
 sum += ic;
 pwm.push_back(new_col);
 }
 }
 }*/

static void generate_targeted_position_weight_matrix(gsl_rng *gr, size_t n,
    size_t m, double target_ic, vector<vector<double> > &pwm) {
  double ic = 0;
  double sum = 0;
  vector<double> new_col(m, 0.0);
  double ave = 0;
  while (fabs(ave - target_ic) > 0.0001) {
    sum = 0;
    pwm.clear();
    while (pwm.size() < n) {
      ic = 0;
      new_col = generate_dirichlet_random_variables(gr, m);
      for (size_t j = 0; j < new_col.size(); ++j)
        if (new_col[j] > 0)
          ic += new_col[j] * abs(log(new_col[j]) / log(2));
      sum += ic;
      pwm.push_back(new_col);
    }
    ave = sum / pwm.size();
  }
}

static size_t roll_die(const Runif &rng, const vector<double> &cumul) {
  return lower_bound(cumul.begin(), cumul.end(), rng.runif(0.0, 1.0))
      - cumul.begin();
}

static void simulate_motif_occurrence(const Runif &rng,
    const vector<vector<double> > &cumul, string &motif_occurrence) {
  motif_occurrence.clear();
  for (size_t i = 0; i < cumul.size(); ++i)
    motif_occurrence += int2base(roll_die(rng, cumul[i]));
}

static void generate_diagnostic_events(const Runif &rng,
    const size_t geo_delta_param, const double geo_rate_param,
    const size_t occurrence_position, const size_t region_size,
    const size_t n_events, vector<size_t> &des) {
  for (size_t i = 0; i < n_events; ++i) {
    const int orig_geo = generate_geometric_random_variable(
        rng, geo_rate_param);
    const int position = occurrence_position + geo_delta_param
        + ((rng.runif(0.0, 1.0) < 0.5) ? orig_geo : -orig_geo);
    des.push_back(
        std::max(0, std::min(position, static_cast<int>(region_size))));
  }

  /* vector<size_t> stars(20,0);
   for (size_t i=0; i<n_events; i+=1)
   stars[des[i] - (int)occurrence_position + 10] +=1;
   for (int i= 0;i<20;i+=1) {
   cerr << i-10 << ":";
   if (stars[i]>0)
   for (size_t j= 0; j< stars[i]; j+=1)
   cerr << "*";
   cerr << endl;
   }*/
}

static GenomicRegion make_diagnostic_event(const GenomicRegion r,
    const size_t de_pos) {
  return GenomicRegion(
      r.get_chrom(), r.get_start() + de_pos, r.get_start() + de_pos + 1,
      r.get_name(), 0, r.get_strand());
}

static void simulate_rnas(const size_t n_des_per_region,
    const size_t geo_delta_param, const double geo_rate_param, const Runif &rng,
    const vector<GenomicRegion> &regions, const vector<vector<double> > cumul,
    vector<GenomicRegion> &diagnostic_events, vector<string> &sequences,
    vector<size_t> &indicators) {

  // const size_t n_current_des = total_des/regions.size();
  const size_t motif_width = cumul.size();

  for (size_t i = 0; i < regions.size(); ++i) {

    const size_t occurrence_position = rng.runif(
        0ul, sequences[i].length() - motif_width);

    string motif_occurrence;
    simulate_motif_occurrence(rng, cumul, motif_occurrence);

    copy(
        motif_occurrence.begin(), motif_occurrence.end(),
        sequences[i].begin() + occurrence_position);

    indicators.push_back(occurrence_position + 1);
//    cerr << occurrence_position + 1 << "\t" << motif_occurrence << endl;

    vector<size_t> current_des;
    generate_diagnostic_events(
        rng, geo_delta_param, geo_rate_param, occurrence_position,
        regions[i].get_width(), n_des_per_region, current_des);

    for (size_t j = 0; j < current_des.size(); ++j)
      diagnostic_events.push_back(
          make_diagnostic_event(regions[i], current_des[j]));
  }
}

static bool check_consistent(const vector<GenomicRegion> &regions,
    const vector<string> &names, const vector<string> &sequences) {
  if (regions.size() != names.size())
    return false;

  for (size_t i = 0; i < regions.size(); ++i)
    if (regions[i].get_name() != names[i]
        || sequences[i].length() != regions[i].get_width())
      return false;
  return true;
}

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

static string pwm2dme(const vector<vector<double> > &M) {

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
      ss << (int) (M[j][b] * 100) << "\t";
    ss << (int) (M[j][alphabet_size - 1] * 100) << endl;
  }
  ss << "XX" << endl;
  ss << "//" << endl;

  return ss.str();
}

static string print_model(const vector<vector<double> > &M, const double p,
    const int delta, const vector<GenomicRegion> &regions,
    const vector<string> &sequences, const vector<size_t> &indicators) {

  // Check that the parameters make sense
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

  ss << "AC\tDE_MODEL" << endl;
  ss << "XX" << endl;
  ss << "TY\tMotif" << endl;
  ss << "XX" << endl;
  ss << "P0\tA\tC\tG\tT" << endl;

  for (size_t j = 0; j < M.size(); j++) {
    ss << "0" << j + 1 << "\t";
    for (size_t b = 0; b < alphabet_size - 1; b++)
      ss << (int) (M[j][b] * N) << "\t";
    ss << (int) (M[j][alphabet_size - 1] * N) << endl;
  }
  ss << "XX" << endl;
  ss << "AT\tGEO_P=" << p << endl;
  ss << "AT\tGEO_delta=" << delta << endl;
  ss << "XX" << endl;

  for (size_t n = 0; n < N; n++) {
    ss << "BS\t" << sequences[n].substr(indicators[n] - 1, M.size()) << "; "
        << regions[n].get_chrom() << ":" << regions[n].get_start() << "-"
        << regions[n].get_end() << "; " << indicators[n] << "; " << M.size()
        << ";  ;" << ((regions[n].get_strand() == '+') ? "p; " : "n;") << endl;
  }
  ss << "XX" << endl;
  ss << "//" << endl;

  return ss.str();
}

int main(int argc, const char **argv) {

  try {

    size_t motif_width = 0;
    size_t n_des_per_region = 10;
    size_t geo_delta_param = 0;
    double geo_rate_param = 0.5;
    int random_number_seed = numeric_limits<int>::max();
    bool VERBOSE = false;
    string pwm_filename = "NA";
    double target_ic = 0;

    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse("simrna", "simulate rna for motif discovery"
        "<pwm> <fasta> <bed> <output>");
    opt_parse.add_opt("width", 'w', "motif width", true, motif_width);
    opt_parse.add_opt(
        "load", 'p', "Name of PWM file (default: random generation)", false,
        pwm_filename);
    opt_parse.add_opt(
        "n-des", 'n', "number of de per region", true, n_des_per_region);
    opt_parse.add_opt(
        "ic", 't', "target average information content", false, target_ic);
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
    if (pwm_filename == "NA" && target_ic == 0) {
      cerr << "For generating PWM, a target information content is needed"
          << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 3) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }

    const string input_fasta_file(leftover_args.front());
    const string input_bed_file(leftover_args[1]);
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<GenomicRegion> regions;
    ReadBEDFile(input_bed_file, regions);

    vector<string> names, sequences;
    read_fasta_file(input_fasta_file, names, sequences);

    if (!check_consistent(regions, names, sequences))
      throw SMITHLABException("inconsistent region names in FASTA/BED files");

    const Runif rng(random_number_seed);
    gsl_rng *gr = gsl_rng_alloc(gsl_rng_mt19937);

    vector<vector<double> > pwm;

    if (pwm_filename != "NA")
      read_position_weight_matrix(pwm_filename, pwm);
    else
      generate_targeted_position_weight_matrix(
          gr, motif_width, smithlab::alphabet_size, target_ic, pwm);

    vector<vector<double> > motif = pwm;

    string dout = pwm2dme(motif);
//    cerr << dout << endl;

    for (size_t i = 0; i < pwm.size(); ++i)
      std::partial_sum(pwm[i].begin(), pwm[i].end(), pwm[i].begin());

    vector<size_t> indicators;
    vector<GenomicRegion> diagnostic_events;

    simulate_rnas(
        n_des_per_region, geo_delta_param, geo_rate_param, rng, regions, pwm,
        diagnostic_events, sequences, indicators);

    // Output the diagnostic events BED file and the updated sequences
    std::ofstream de_out(string(outfile + ".bed").c_str());
    copy(
        diagnostic_events.begin(), diagnostic_events.end(),
        std::ostream_iterator<GenomicRegion>(de_out, "\n"));

    std::ofstream seqs_out(string(outfile + ".fa").c_str());
    for (size_t i = 0; i < sequences.size(); ++i)
      seqs_out << '>' << regions[i].get_name() << endl << sequences[i] << endl;

    cout
        << print_model(
            motif, geo_rate_param, geo_delta_param, regions, sequences,
            indicators);

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
