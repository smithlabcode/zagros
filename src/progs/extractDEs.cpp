/**
  \file extractDEs.cpp
  \brief TODO

  \authors Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

  \section copyright Copyright Details
  Copyright (C) 2014
  University of Southern California,
  Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

  \section license License Details
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

  \section bugs Known Bugs

  \section history Revision History
**/

#include <string>
#include <vector>
#include <tr1/unordered_map>

#include "IO.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::ostream;
using std::endl;
using std::pair;
using std::tr1::unordered_map;
using std::stringstream;

//


/**
 * \brief check whether a set of genomic intervals is free of overlaps. The
 *        intervals need not be sorted, but they must all be from the same
 *        chromosome (an exception is thrown if they aren't). Since genomic
 *        regions don't include their end coordinate, we don't consider two
 *        genomic intervals to overlap if one of them has a start coordinate
 *        equal to the end coordinate of the other.
 * \param   target_chrom  the set of genomic regions to test
 * \return  true if the set is free of overlaps, false otherwise
 * \throw   SMITHLABException if the regions don't all come from the
 *          same chromosome
 * \todo    there are now several copies of slightly different code in zagros
 *          for checking whether a set of genomic regions has any overlapping
 *          elements; they need to be consolidated.
 */
static bool
check_overlapping_chrom(const vector<GenomicRegion> &targets_chrom) {
  bool first = true;
  string chrom = "";
  vector<pair<size_t, bool> > boundaries;
  for (size_t i = 0; i < targets_chrom.size(); ++i) {
    if (first) {
      first = false;
      chrom = targets_chrom[i].get_chrom();
    } else {
      if (chrom != targets_chrom[i].get_chrom()) {
        stringstream ss;
        ss << "failed checking overlap of chromosomes. Supplied regions are "
           << "on multiple chromosomes: " << chrom << " and "
           << targets_chrom[i].get_chrom();
        throw SMITHLABException(ss.str());
      }
    }
    boundaries.push_back(std::make_pair(targets_chrom[i].get_start(), false));
    boundaries.push_back(std::make_pair(targets_chrom[i].get_end()-1, true));
  }
  sort(boundaries.begin(), boundaries.end());

  size_t count = 0;
  for (size_t i = 0; i < boundaries.size(); ++i)
    if (boundaries[i].second)
      --count;
    else {
      ++count;
      if (count > 1) {
        cerr << targets_chrom[i].tostring() << endl;
        return true;
      }
    }
  return false;
}

static bool
check_overlapping(const vector<GenomicRegion> &targets) {

  vector<vector<GenomicRegion> > targets_chroms;
  separate_chromosomes(targets, targets_chroms);
  bool ret_val = false;
  for (size_t i = 0; i < targets_chroms.size() && !ret_val; ++i)
    ret_val |= check_overlapping_chrom(targets_chroms[i]);
  return ret_val;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////
/////////////  STUFF FOR DIAGNOSTIC EVENTS
/////////////

struct DE {
  std::string type;
  size_t position;
  std::string base;
  std::string read;
};

static size_t
convertString(const string& s) {
  std::istringstream buffer(s);
  size_t value;
  buffer >> value;
  return value;
}


static void
sift_single_chrom(vector<GenomicRegion> &other_regions,
                  vector<GenomicRegion> &regions,
                  vector<GenomicRegion> &good_regions) {

  typedef vector<GenomicRegion>::iterator region_itr;

  for (size_t i = 0; i < regions.size(); ++i) {
    region_itr closest(find_closest(other_regions, regions[i]));
    if (closest->overlaps(regions[i])) {
      good_regions.push_back(regions[i]);
      good_regions.back().set_name(closest->get_name());
    }
  }
}

static void
sift(vector<GenomicRegion> &other_regions,
     vector<GenomicRegion> &regions) {

  vector<vector<GenomicRegion> > other_regions_by_chrom;
  separate_chromosomes(other_regions, other_regions_by_chrom);

  unordered_map<string, size_t> chrom_lookup;
  for (size_t i = 0; i < other_regions_by_chrom.size(); ++i)
    chrom_lookup[other_regions_by_chrom[i].front().get_chrom()] = i;
  const vector<GenomicRegion> dummy;

  vector<vector<GenomicRegion> > regions_by_chrom;
  separate_chromosomes(regions, regions_by_chrom);
  regions.clear();

  vector<GenomicRegion> good_regions;
  for (size_t i = 0; i < regions_by_chrom.size(); ++i) {
    const unordered_map<string, size_t>::const_iterator j =
      chrom_lookup.find(regions_by_chrom[i].front().get_chrom());
    if (j != chrom_lookup.end()) {
      sift_single_chrom(other_regions_by_chrom[j->second],
            regions_by_chrom[i], good_regions);
    }

  }
  regions.swap(good_regions);

  unordered_map<string, size_t> name_lookup;
  for (size_t i = 0; i < other_regions.size(); ++i)
    name_lookup[other_regions[i].get_name()] = i;

  string prev_name = "";
  for (size_t i = 0; i < regions.size(); ++i)
    if (regions[i].get_name() != prev_name) {
      const size_t idx(name_lookup[regions[i].get_name()]);
      other_regions[idx].set_strand(regions[i].get_strand());
      prev_name = regions[i].get_name();
    }
}

static void
parse_diagnostic_events(const string &de_string,
                        vector<DE> &des) {
  vector<string> parts(smithlab::split_whitespace_quoted(de_string));
  if (parts.size() > 0)
    for (size_t i = 0; i < parts.size(); ++i) {
      DE de;
      size_t index = parts[i].find(">");
      if (index != string::npos) {
        de.type = "mutation";
        de.position = convertString(parts[i].substr(0, index - 1));
        de.base = parts[i].substr(index - 1, 1);
        de.read = parts[i].substr(index + 1, 1);
      }
      else {
        index = parts[i].find("+");
        if (index != string::npos) {
          de.type = "deletion";
          de.position = convertString(parts[i].substr(0, index));
          de.base = parts[i].substr(index + 1, parts[i].length() - index - 1);
        }
    else {
          index = parts[i].find("-");
          if (index != string::npos) {
            de.type = "insertion";
            de.position = convertString(parts[i].substr(0, index));
            de.base = parts[i].substr(index + 1, parts[i].length() - index - 1);
          }
        }
      }
      des.push_back(de);
    }
}

static void
add_diagnostic_events_iCLIP(const MappedRead &region,
                            vector<GenomicRegion> &diagnostic_events) {

  vector<DE> des;
  const string name("X");
  parse_diagnostic_events(region.scr, des);
  bool contains_de = true;
  if (des.size() > 0)
    for (size_t i = 0; i < des.size(); ++i)
      if ((region.r.pos_strand() && des[i].type == "insertion"
       && des[i].base == "T")
          || (!region.r.pos_strand() && des[i].type == "insertion"
              && des[i].base == "A"))
        contains_de = false;
  if (contains_de) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events.back().get_name() == name)
        score = diagnostic_events.back().get_score() + 1;
    if (region.r.pos_strand()) {
      GenomicRegion gr(region.r.get_chrom(),
               region.r.get_start(),
               region.r.get_start() + 1,
               name, score, region.r.get_strand());
      diagnostic_events.push_back(gr);
    }
    else {
      GenomicRegion gr(region.r.get_chrom(), region.r.get_end(),
               region.r.get_end() + 1,
               name, score, region.r.get_strand());
      diagnostic_events.push_back(gr);
    }
  }
}

static void
add_diagnostic_events_hCLIP(const MappedRead &region,
                            vector<GenomicRegion> &diagnostic_events) {

  vector<DE> des;
  const string name("X");
  parse_diagnostic_events(region.scr, des);
  if (des.size() == 1) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events.back().get_name() == name)
        score = diagnostic_events.back().get_score() + 1;
    for (size_t j = 0; j < des.size(); ++j) {
      if (region.r.pos_strand() &&
      des[j].type == "insertion" && des[j].base == "T") {
        GenomicRegion gr(region.r.get_chrom(),
             region.r.get_start() + des[j].position - 1,
             region.r.get_start() + des[j].position, name, score,
             region.r.get_strand());
        diagnostic_events.push_back(gr);
      }
      else if (!region.r.pos_strand() &&
           des[j].type == "insertion" && des[j].base == "A") {
        GenomicRegion gr(region.r.get_chrom(),
             region.r.get_start() + des[j].position - 2,
             region.r.get_start() + des[j].position - 1, name, score,
             region.r.get_strand());
        diagnostic_events.push_back(gr);
      }
    }
  }
}

static void
add_diagnostic_events_pCLIP(const MappedRead &region,
                            vector<GenomicRegion> &diagnostic_events) {
  vector<DE> des;
  parse_diagnostic_events(region.scr, des);
  const string name("X");
  if (des.size() == 1) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events.back().get_name() == name)
        score = diagnostic_events.back().get_score() + 1;
    for (size_t j = 0; j < des.size(); ++j) {
      if (region.r.pos_strand() && des[j].type == "mutation"
          && des[j].base == "T" && des[j].read == "C") {
        GenomicRegion gr(region.r.get_chrom(),
             region.r.get_start() + des[j].position - 1,
             region.r.get_start() + des[j].position, name, score,
             region.r.get_strand());
        diagnostic_events.push_back(gr);
      } else if (!region.r.pos_strand() && des[j].type == "mutation"
         && des[j].base == "A" && des[j].read == "G") {
        GenomicRegion gr(region.r.get_chrom(),
             region.r.get_start() + des[j].position - 2,
             region.r.get_start() + des[j].position - 1, name, score,
             region.r.get_strand());
        diagnostic_events.push_back(gr);
      }
    }
  }
}

static void
add_diagnostic_events(const MappedRead &region,
                      const string experiment,
                      vector<GenomicRegion> &diagnostic_events) {
  if (experiment == "iCLIP")
    add_diagnostic_events_iCLIP(region, diagnostic_events);
  else if (experiment == "hCLIP")
    add_diagnostic_events_hCLIP(region, diagnostic_events);
  else if (experiment == "pCLIP")
    add_diagnostic_events_pCLIP(region, diagnostic_events);
  else
    throw SMITHLABException("The technology was not recognized!");
}

static void
add_diagnostic_events(const vector<MappedRead> &regions,
                      const string experiment,
                      vector<GenomicRegion> &diagnostic_events) {

  for (size_t i = 0; i < regions.size(); ++i) {
    add_diagnostic_events(regions[i], experiment, diagnostic_events);
  }
}


static void
separate_mapped_reads_chromosomes(const vector<MappedRead> &regions,
                                  vector<vector<MappedRead> > &separated_by_chrom) {
  typedef unordered_map<string, vector<MappedRead> > Separator;
  Separator separator;
  for (vector<MappedRead>::const_iterator i = regions.begin();
       i != regions.end(); ++i) {
    const string the_chrom(i->r.get_chrom());
    if (separator.find(the_chrom) == separator.end())
      separator[the_chrom] = vector<MappedRead>();
    separator[the_chrom].push_back(*i);
  }
  separated_by_chrom.clear();
  for (Separator::iterator i = separator.begin(); i != separator.end(); ++i)
    separated_by_chrom.push_back(i->second);
}

static void
make_inputs(const vector<MappedRead> &mapped_reads,
            const string experiment,
            vector<GenomicRegion> &diagnostic_events) {
  vector<vector<MappedRead> > separated_by_chrom;
  separate_mapped_reads_chromosomes(mapped_reads, separated_by_chrom);
  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    add_diagnostic_events(separated_by_chrom[i], experiment, diagnostic_events);
    separated_by_chrom[i].clear();
  }
}


/***
 * \summary Count the number of DEs that fall at each position within a
 *          set of genomic regions
 * \param events    must be sorted
 * \param regions   must be non-overlapping and sorted.
 * \param binCounts \todo
 * \throws SMITHLABException \todo
 */
void
countDEs (const vector<GenomicRegion> &events,
          const vector<GenomicRegion> &regions, vector<vector<size_t> > &counts) {
  if (!check_sorted(events))
    throw SMITHLABException("DE event locations must be in sorted order");
  if (!check_sorted(regions))
    throw SMITHLABException("Regions must be in sorted order");
  if (check_overlapping(regions))
    throw SMITHLABException("Regions must be non-overlapping");

  if ((regions.size() == 0) || (events.size() == 0)) return;

  counts.clear();
  counts.resize(regions.size());
  for (size_t i = 0; i < regions.size(); ++i)
    counts[i].resize(regions[i].get_width());

  size_t bin_idx = 0, de_idx = 0;
  while (bin_idx < regions.size()) {
    // while we haven't run out of des, and the de chrom is less than the bin,
    // or it's the same but the end of the de is still less than the start of
    // the bin...
    while ((de_idx < events.size()) &&
           ((events[de_idx].get_chrom() < regions[bin_idx].get_chrom()) ||
            ((events[de_idx].get_chrom() == regions[bin_idx].get_chrom()) &&
                (events[de_idx].get_end() < regions[bin_idx].get_start())))) {
      de_idx += 1;
    }
    // if not (we ran out of des, or the de chrom is greater than the bin, or
    // it's the same but the start of the de is greater than the end of the
    // bin...)
    if (!((de_idx == events.size()) ||
        (events[de_idx].get_chrom() > regions[bin_idx].get_chrom()) ||
        ((events[de_idx].get_chrom() == regions[bin_idx].get_chrom()) &&
            (events[de_idx].get_start() > regions[bin_idx].get_end())))) {
      // the de falls inside the bin; count them up while they
      // continue to do so
      while ((de_idx < events.size()) &&
             (regions[bin_idx].overlaps(events[de_idx]))) {
        size_t relIndex = events[de_idx].get_start() -\
                          regions[bin_idx].get_start();
        assert (relIndex < counts[bin_idx].size());
        counts[bin_idx][relIndex] += 1;
        de_idx += 1;
      }
    }
    bin_idx += 1;
  }
}


/****
 * \summary TODO
 */
int
main(int argc, const char **argv) {
  try {
    string outputFn = "";
    string regionsFn = "";
    string mapper = "rmap";
    string tech = "iCLIP";
    bool VERBOSE = false;


    // TODO specify defaults below for options

    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse(strip_path(argv[0]),
                                  "program for extracting diagnostic events "
                                  "from mapped reads files");
    opt_parse.add_opt("output", 'o', "Write output to this file (STDOUT if "
                      "omitted).", OptionParser::OPTIONAL, outputFn);
    opt_parse.add_opt("regions", 'r', "the genomic regions, in BED format, "
                      "corresponding to the input sequences for Zagros.",
                      OptionParser::REQUIRED, regionsFn);
    opt_parse.add_opt("mapper", 'm', "the mapper used to map the reads",
                       OptionParser::OPTIONAL, mapper);
    opt_parse.add_opt("tech", 't', "the technology type used in the experiment",
                      OptionParser::OPTIONAL, tech);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      OptionParser::OPTIONAL, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    if (leftover_args.empty()) {
      throw SMITHLABException("Must provide mapped reads file");
    }
    if (leftover_args.size() > 1) {
      stringstream ss;
      ss << "Multiple mapped reads files detected: ";
      for (size_t i = 0; i < leftover_args.size(); i++) {
        ss << leftover_args[i];
        if (i != leftover_args.size() - 1) ss << ", ";
        else ss << "; can only handle one at a time";
      }

      throw SMITHLABException(ss.str());
    }
    const string mappedReadFn(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE) cerr << "LOADING MAPPING INFORMATION... ";
    vector<MappedRead> mapped_reads;
    load_mapped_reads(mappedReadFn, mapper, mapped_reads);
    if (VERBOSE) cerr << "DONE" << endl;

    if (VERBOSE) cerr << "LOADING TARGET REGIONS... ";
    vector<GenomicRegion> targets;
    ReadBEDFile(regionsFn, targets);
    if (check_overlapping(targets))
      throw SMITHLABException("Target regions are overlapping!");
    if (VERBOSE) cerr << "DONE" << endl;

    if (VERBOSE) cerr << "PROCESSING MAPPED READS... ";
    vector<GenomicRegion> de_regions;
    make_inputs(mapped_reads, tech, de_regions);
    sort(de_regions.begin(), de_regions.end());
    sift(targets, de_regions);
    sort(de_regions.begin(), de_regions.end());
    vector<vector<size_t> > counts;
    countDEs(de_regions, targets, counts);
    if (VERBOSE) cerr << "DONE" << endl;

    if (VERBOSE) cerr << "WRITING OUTPUT... ";
    std::ofstream of;
    if (!outputFn.empty()) of.open(outputFn.c_str());
    ostream ostrm(outputFn.empty() ? cout.rdbuf() : of.rdbuf());
    if (!ostrm.good()) {
      if (!outputFn.empty())
        throw SMITHLABException("Failed to open " + outputFn + " for writing");
      else
        throw SMITHLABException("Unable to write to standard out");
    }
    for (size_t i = 0; i < counts.size(); ++i) {
      for (size_t j = 0; j < counts[i].size(); ++j) {
        ostrm << counts[i][j];
        if (j != counts[i].size() - 1) ostrm << ", ";
      }
      ostrm << endl;
    }
    if (VERBOSE) cerr << endl;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}









