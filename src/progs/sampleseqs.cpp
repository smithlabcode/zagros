#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "RNG.hpp"
#include "OptionParser.hpp"

using std::numeric_limits;
using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;

static void sample_seqs(const Runif &rng, const vector<GenomicRegion> &regions,
    const vector<string> &sequences, vector<GenomicRegion> &sample_regions,
    vector<string> &sample_sequences, const size_t n) {

  // const size_t n_current_des = total_des/regions.size();
  //
  for (size_t i = 0; i < n; ++i) {
    const size_t sequence_number = rng.runif(0ul, sequences.size());
    sample_sequences.push_back(sequences[sequence_number]);
    sample_regions.push_back(regions[sequence_number]);
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

int main(int argc, const char **argv) {

  try {

    size_t n_seqs = 0;
    int random_number_seed = numeric_limits<int>::max();
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse("sample_sequence", "sample n rna sequences"
        "<fasta> <bed> <output>");
    opt_parse.add_opt("n-seqs", 'n', "number of sequences", true, n_seqs);
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

    vector<GenomicRegion> sample_regions;
    vector<string> sample_sequences;
    sample_seqs(
        rng, regions, sequences, sample_regions, sample_sequences, n_seqs);

    // Output the diagnostic events BED file and the updated sequences
    std::ofstream regions_out(string(outfile + ".bed").c_str());
    copy(
        sample_regions.begin(), sample_regions.end(),
        std::ostream_iterator<GenomicRegion>(regions_out, "\n"));

    std::ofstream seqs_out(string(outfile + ".fa").c_str());
    for (size_t i = 0; i < sample_sequences.size(); ++i)
      seqs_out << '>' << sample_regions[i].get_name() << endl
          << sample_sequences[i] << endl;

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
