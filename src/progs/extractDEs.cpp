/**
  \file extractDEs.cpp
  \brief This is a program that is part of the Zagros package. It's function
         is to compute the counts of diagnostic events within a set of target
         regions (that the user specifies) by extracting them from mapping
         results.

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

// STL includes
#include <string>
#include <vector>
#include <tr1/unordered_map>

// Zagros common includes
#include "IO.hpp"
#include "IntervalTree.hpp"

// Smithlab common icnludes
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::vector;
using std::ostream;
using std::stringstream;
using std::tr1::unordered_map;


/******************************************************************************
 *      structures, and functions for manipulating diagnostic events
 *****************************************************************************/

/**
 * \brief The DERoi structure defines a 'region of interest' to extract counts
 *        of diagnostic events in.
 */
struct DERoi {
public:
  DERoi(GenomicRegion &r_p, size_t rank_p) : r(r_p), rank(rank_p) {
    this->counts.resize(r.get_width());
  }
  GenomicRegion r;          /** the region covered **/
  vector<size_t> counts;    /** counts[i] is the num. of DEs at rel. pos i **/
  size_t rank;              /** position in input; used for sorting **/
};

/**
 * \brief function for getting the start coordinate of a DERoi from a pointer;
 *        this is just used for the IntervalTree, which takes a function
 *        pointer for how to extract a start coordinate from an object of
 *        type T (here pointer to DERoi)
 */
size_t getDERoiStart(DERoi *a) {return a->r.get_start();}

/**
 * \brief function for getting the end coordinate of a DERoi from a pointer;
 *        this is just used for the IntervalTree, which takes a function
 *        pointer for how to extract an end coordinate from an object of
 *        type T (here pointer to DERoi)
 */
size_t getDERoiEnd(DERoi *a) {return a->r.get_end();}

/**
 * \brief a functor used for sorting DERoi objects; just imposes an ordering
 *        base on the rank (position within input file) or the ROI.
 */
struct DERoiRankComparator {
  inline bool operator()(const DERoi *a, const DERoi *b) {
    return a->rank < b->rank;
  }
};

/**
 * \brief This structure wraps up information about a single diagnostic event.
 */
struct DE {
  std::string type;         /** type of DE. insertion/deletion/mutation   **/
  size_t position;          /** location of DE relative to read start pos **/
  std::string refBase;      /** the base in the reference (genome)        **/
  std::string readBase;     /** the base in the read                      **/
};

/******************************************************************************
 *            simple helper functions for manipulating types
 *****************************************************************************/

/**
 * \brief convert a string to a size_t
 * \param   s   the string to convert
 * \return      the size_t corresponding to s
 */
static size_t
stringToSize_t(const string& s) {
  std::istringstream buffer(s);
  size_t value;
  buffer >> value;
  return value;
}

/**
 * \brief convert a bool to a string representation
 * \param b the bool to convert
 * \return the string "true" if b is true, otherwise the string "false"
 */
static string
toString(const bool b) {
  if (b) return "true";
  return "false";
}

/******************************************************************************
 *   functions for doing the actual extraction of DEs from mapping results
 *****************************************************************************/


/**
 * \brief parse a Bowtie mismatch string to extract any diagnostic events
 * 
 *        Bowtie (native format) mapped reads have a tab-delimited
 *        format; the eighth field is the "mismatch string". Note that only 
 *        this field should be provided to this function, not the full read string. 
 *        The mismatch string is comma delimited; we're interested in entries 
 *        that specify insertions, deletions or mismatches between the read and 
 *        the reference (there may be other entries, but we don't care about them).
 *        Mismatch:  Format is 'offset':'refbase'>'readbase'
 *        Insertion: Format is 'offset'+'insertedbases'
 *        Deletion:  Format is 'offset'-'refbase'
 *        --
 *        deletion means the base(s) is (are) in the ref., not in the read
 *        insertion means the base(s) is (are) in the read, not the ref.
 *        offset is the mm/indel offset from the mapping location (not the
 *        position within the read).
 *
 * \param bowtie_mmString   the mismatch string from a Bowtie read
 *                          (can be the empty string, or all whitespace)
 * \param des               any extracted DEs will be added to the back of this
 *                          vector. No existing entries in the vector will be
 *                          modified.
 **/
static void
extractDEs_bowtie(const string &bowtie_mmString, vector<DE> &des) {
  vector<string> parts(smithlab::split_whitespace_quoted(bowtie_mmString));
  for (size_t i = 0; i < parts.size(); ++i) {
    DE de;
    size_t idx = parts[i].find(">");
    if (idx != string::npos) {
      de.type = "mutation";
      de.position = stringToSize_t(parts[i].substr(0, idx - 1));
      de.refBase = parts[i].substr(idx - 1, 1);
      de.readBase = parts[i].substr(idx + 1, 1);
      des.push_back(de);
    }
    else {
      idx = parts[i].find("+");
      if (idx != string::npos) {
        de.type = "insertion";
        de.position = stringToSize_t(parts[i].substr(0, idx));
        de.readBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
        de.refBase = "-";
        des.push_back(de);
      } else {
        idx = parts[i].find("-");
        if (idx != string::npos) {
          de.type = "deletion";
          de.position = stringToSize_t(parts[i].substr(0, idx));
          de.refBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
          de.readBase = "-";
          des.push_back(de);
        }
      }
    }
  }
}


/**
 * \brief parse a Novoalign mismatch string to extract any diagnostic events
 *
 *        Novoalign (native format) mapped reads have a tab-delimited
 *        format; there are a lot of possible fields, and they vary depending
 *        on the exact input provided, but the final one is what I call here
 *        the "mismatch string". Note that only this field should be provided
 *        to this function, not the full read string. The mismatch string is
 *        space delimited; we're interested in entries that specify insertions,
 *        deletions or mismatches between the read and the reference (there may
 *        be other entries, but we don't care about them).
 *        Mismatch:  Format is 'offset''refbase'>'readbase'
 *        Insertion: Format is 'offset'+'insertedbases'
 *        Deletion:  Format is 'offset'-'refbase'
 *        --
 *        deletion means the base(s) is (are) in the ref., not in the read
 *        insertion means the base(s) is (are) in the read, not the ref.
 *        offset is the mm/indel offset from the mapping location (not the
 *        position within the read).
 *
 * \param novo_mmString     the mismatch string from a Novoalign read
 *                          (can be the empty string, or all whitespace)
 * \param des               any extracted DEs will be added to the back of this
 *                          vector. No existing entries in the vector will be
 *                          modified.
 */
static void
extractDEs_novo(const string &novo_mmString, vector<DE> &des) {
  vector<string> parts(smithlab::split_whitespace_quoted(novo_mmString));
  for (size_t i = 0; i < parts.size(); ++i) {
    DE de;
    size_t idx = parts[i].find(">");
    if (idx != string::npos) {
      de.type = "mutation";
      de.position = stringToSize_t(parts[i].substr(0, idx - 1));
      de.refBase = parts[i].substr(idx - 1, 1);
      de.readBase = parts[i].substr(idx + 1, 1);
      des.push_back(de);
    }
    else {
      idx = parts[i].find("+");
      if (idx != string::npos) {
        de.type = "insertion";
        de.position = stringToSize_t(parts[i].substr(0, idx));
        de.readBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
        de.refBase = "-";
        des.push_back(de);
      } else {
        idx = parts[i].find("-");
        if (idx != string::npos) {
          de.type = "deletion";
          de.position = stringToSize_t(parts[i].substr(0, idx));
          de.refBase = parts[i].substr(idx + 1, parts[i].length() - idx - 1);
          de.readBase = "-";
          des.push_back(de);
        }
      }
    }
  }
}

/**
 * \brief given a mapped read from an iCLIP experiment, produced by a
 *        particular type of short read mapper, extract the diagnostic events
 *        and add their genomic loci to the provided vector. An iCLIP read is
 *        considered to have a diagnostic event as long as there are no 'T'
 *        deletions in the read. A 'T' deletion would imply that the reverse
 *        transcriptase read through the cross-link location, producing the 'T'
 *        deletion, and failed to terminate at it. Any iCLIP read satisfying
 *        that requirement is considered to have a DE at its 5' end.
 * \param read the mapped read
 * \param mapper        The name of the short read mapper used
 * \param diagEventLocs the genomic loci of the extracted diagnostic events
 *                      will be added to the end of this vector. Existing
 *                      elements in the vector are left unchanged.
 */
static void
addDiagEvents_iCLIP(const MappedRead &read, const string &mapper,
                    vector<GenomicRegion> &diagEventLocs) {
  vector<DE> des;
  if (mapper == "novoalign") 
    extractDEs_novo(read.scr, des);
  else if (mapper == "bowtie")
    extractDEs_bowtie(read.scr, des);
  else if (mapper != "rmap")
    throw SMITHLABException("unsupported short-read mapper for iCLIP "
                               "data: " + mapper);

  // check for 'T' deletion in the read that signifies an RT read-through of
  // the cross-link location; we consider such reads to not contain DEs.
  bool contains_de = true;
  for (size_t i = 0; i < des.size(); ++i) {
    const bool posStrnd = read.r.pos_strand();
    if ((posStrnd && des[i].type == "deletion" && des[i].refBase == "T") ||
        (!posStrnd && des[i].type == "deletion" && des[i].refBase == "A")) {
      contains_de = false;
    }
  }

  if (contains_de) {
    const string name("X");
    const size_t score = 0;
    if (read.r.pos_strand()) {
      GenomicRegion g(read.r.get_chrom(), read.r.get_start(),
                      read.r.get_start() + 1, name, score, read.r.get_strand());
      diagEventLocs.push_back(g);
    }
    else {
      GenomicRegion g(read.r.get_chrom(), read.r.get_end() - 1,
                       read.r.get_end(), name, score, read.r.get_strand());
      diagEventLocs.push_back(g);
    }
  }
}


/**
 * \brief given a mapped read from a HITS-CLIP experiment, produced by a
 *        particular type of short read mapper, extract any diagnostic events
 *        from it and add their genomic loci to the  provided vector.
 *        A diagnostic event in a HITS-CLIP read is the deletion of a T. We
 *        only consider reads with one DE; x DEs implies x-1 are spurious, and
 *        we wouldn't know which was the correct one, so we discard such cases.
 * \param read          The mapped read
 * \param mapper        The name of the short read mapper used
 * \param diagEventLocs the genomic loci of the extracted diagnostic events
 *                      will be added to the end of this vector. Existing
 *                      elements in the vector are left unchanged.
 */
static void
addDiagEvents_hCLIP(const MappedRead &read, const string &mapper,
                    vector<GenomicRegion> &diagEventLocs) {
  vector<DE> des;
  if (mapper == "novoalign")
    extractDEs_novo(read.scr, des);
  else if (mapper == "bowtie")
    extractDEs_bowtie(read.scr, des);
  else
    throw SMITHLABException("unsupported short-read mapper for hCLIP "
                               "data: " + mapper);

//  extractDEs_novo(read.scr, des);
  if (des.size() == 1) {
    const string name("X");
    const size_t score = 0;
    if (read.r.pos_strand() &&
    des[0].type == "deletion" && des[0].refBase == "T") {
      GenomicRegion gr(read.r.get_chrom(),
           read.r.get_start() + des[0].position - 1,
           read.r.get_start() + des[0].position, name, score,
           read.r.get_strand());
      diagEventLocs.push_back(gr);
//      cout << gr.get_start() << endl;
    }
    else if (!read.r.pos_strand() &&
         des[0].type == "deletion" && des[0].refBase == "A") {
      GenomicRegion gr(read.r.get_chrom(),
           read.r.get_start() + des[0].position - 2,
           read.r.get_start() + des[0].position - 1, name, score,
           read.r.get_strand());
      diagEventLocs.push_back(gr);
    }
  }
}

/**
 * \brief given a read from a PAR-CLIP experiment, produced by a particular
 *        type of short-read mapper, extract any diagnostic events from it
 *        and add their genomic loci to the provided vector. A DE for PAR-CLIP
 *        is a T->C mutation in the read (i.e. read has 'C', ref. genome has
 *        'T'). We only consider reads with one DE; x DEs implies x-1 are
 *        spurious, and we wouldn't know which was the correct one, so we
 *        discard such cases.
 * \param read          The mapped read
 * \param mapper        The name of the short read mapper used
 * \param diagEventLocs the genomic loci of the extracted diagnostic events
 *                      will be added to the end of this vector. Existing
 *                      elements in the vector are left unchanged.
 */
static void
addDiagEvents_pCLIP(const MappedRead &read, const string &mapper,
                    vector<GenomicRegion> &diagEventLocs) {
  vector<DE> des;
  if (mapper == "novoalign")
    extractDEs_novo(read.scr, des);
  else if (mapper == "bowtie")
    extractDEs_bowtie(read.scr, des);
  else
    throw SMITHLABException("unsupported short-read mapper for pCLIP "
                               "data: " + mapper);

  if (des.size() == 1) {
    const string name("X");
    const size_t score = 1;
    if (read.r.pos_strand() && des[0].type == "mutation"
        && des[0].refBase == "T" && des[0].readBase == "C") {
      GenomicRegion gr(read.r.get_chrom(),
          read.r.get_start() + des[0].position - 1,
          read.r.get_start() + des[0].position, name, score,
          read.r.get_strand());
      diagEventLocs.push_back(gr);
    } else if (!read.r.pos_strand() && des[0].type == "mutation"
       && des[0].refBase == "A" && des[0].readBase == "G") {
      GenomicRegion gr(read.r.get_chrom(),
          read.r.get_start() + des[0].position - 2,
          read.r.get_start() + des[0].position - 1, name, score,
          read.r.get_strand());
      diagEventLocs.push_back(gr);
    }
  }
}

/**
 * \brief for a given mapped read, for a given experiment type (iCLIP,
 *        PAR-CLIP, HITS-CLIP), from a given short read mapper, extract the
 *        diagnostic events from the read and add their genomic loci to the
 *        provided vector.
 * \param read          extract DEs from these mapped reads
 * \param experiment    Which type of experiment (iCLIP, PAR-CLIP,
 *                      HITS-CLIP)
 * \param mapper        Which short read mapper was used
 * \param diagEvents    the genomic loci of diagnostic events extracted from
 *                      the reads are added to this vector. Any existing
 *                      entries are NOT modified.
 */
static void
addDiagEvents(const MappedRead &read, const string &experiment,
              const string &mapper, vector<GenomicRegion> &diagEvents) {
  if (experiment == "iCLIP") addDiagEvents_iCLIP(read, mapper, diagEvents);
  else if (experiment == "hCLIP") addDiagEvents_hCLIP(read, mapper, diagEvents);
  else if (experiment == "pCLIP") addDiagEvents_pCLIP(read, mapper, diagEvents);
  else throw SMITHLABException("The technology '" + experiment +\
                               "' was not recognized!");
}

/**
 * \brief for a set of mapped reads from a given experiment type (iCLIP,
 *        PAR-CLIP, HITS-CLIP), mapped with a given short read mapper, extract
 *        the diagnostic events from the reads and add their genomic loci to
 *        the provided vector.
 * \param reads          extract DEs from these mapped reads
 * \param experiment     Which type of experiment (iCLIP, PAR-CLIP,
 *                       HITS-CLIP)
 * \param mapper         Which short read mapper was used
 * \param diagEvents     the genomic loci of diagnostic events extracted from
 *                       the reads are added to this vector. Any existing
 *                       entries are NOT modified.
 */
static void
addDiagEvents(const vector<MappedRead> &reads, const string &experiment,
              const string &mapper, vector<GenomicRegion> &diagEvents) {
  for (size_t i = 0; i < reads.size(); ++i)
    addDiagEvents(reads[i], experiment, mapper, diagEvents);
}


/******************************************************************************
 * functions for loading and manipulating interval trees of the target regions
 *****************************************************************************/
typedef unordered_map<string, IntervalTree<DERoi*, size_t> > TreeMap;

/**
 * \summary load a map of interval trees for diag. event regions of interest.
 *          Note that this function allocates its memory on the heap, and the
 *          caller is responsible for deallocating this memory via a
 *          corresponding call to deleteIntervalTrees_DERoi
 * \param regions_fn load regions from this file
 * \param trees      after the call, trees[chr] will contain the interval tree
 *                   for the chromosome chr. All existing entries are cleared.
 */
static void
loadIntervalTrees_DERoi(const string &regions_fn, TreeMap &trees) {
  typedef unordered_map<string, vector<DERoi*> > regionMap;
  trees.clear();
  vector<GenomicRegion> targets;
  ReadBEDFile(regions_fn, targets);
  regionMap byChrom;
  for (size_t i = 0; i < targets.size(); ++i) {
    if (byChrom.count(targets[i].get_chrom()) == 0)
      byChrom[targets[i].get_chrom()] = vector<DERoi*>();
    byChrom[targets[i].get_chrom()].push_back(new DERoi(targets[i],i));
  }
  for (regionMap::iterator it = byChrom.begin(); it != byChrom.end(); ++it) {
    trees.insert(std::make_pair<std::string, IntervalTree<DERoi*, size_t> >
                (it->first, IntervalTree<DERoi*, size_t> (it->second,
                                                 getDERoiStart, getDERoiEnd)));
  }
}

/**
 * \brief clean up a TreeMap by deallocating the memory it reserved on the heap
 *        and removing all of the keys from the map.
 * \param t the TreeMap to delete.
 */
static void
deleteIntervalTrees_DERoi(TreeMap &t) {
  for (TreeMap::iterator it = t.begin(); it != t.end(); it++) {
    vector<DERoi*> v;
    it->second.squash(v);
    for (size_t i = 0; i < v.size(); ++i) delete v[i];
  }
  t.clear();
}

/******************************************************************************
 *               functions for preparing and writing output
 *****************************************************************************/

/**
 * \brief write diagnostic event counts to file.
 * \param counts
 * \param fn                    file to write output to. If this is the empty
 *                              string, we'll write to stdout instead.
 * \param VERBOSE               if true, output additional status messages to
 *                              stderr about the progress of the computation
 * \throw SMITHLABException     if the specified file cannot be opened for
 *                              writing, or if the filename is equal to the
 *                              empty string and we fail to get a good handle
 *                              for stdout.
 */
void
writeDEs(const vector< vector <size_t> > &counts, const string &fn,
         const bool VERBOSE=false) {
  if (VERBOSE) cerr << "WRITING OUTPUT... ";
  std::ofstream of;
  if (!fn.empty()) of.open(fn.c_str());
  ostream ostrm(fn.empty() ? cout.rdbuf() : of.rdbuf());
  if (!ostrm.good()) {
    if (!fn.empty())
      throw SMITHLABException("Failed to open " + fn + " for writing");
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
  if (VERBOSE) cerr << "DONE" << endl;
}

/******************************************************************************
 *              Main entry point and code for extracting DEs
 *****************************************************************************/

/**
 * \brief given a file of mapped reads for a particular type of experiment
 *        from a particular short read mapper, and a file with a set of
 *        target regions in BED format, count the number of diagnostic
 *        events that occur at each position within each of the target
 *        regions.
 * \param mappedReads_fn        the file with the mapped reads
 * \param targetRegions_fn      the file with the target regions (BED format)
 * \param experimentType        Which type of experiment (iCLIP, PAR-CLIP,
 *                              HITS-CLIP)
 * \param mapperType            Which short read mapper was used
 * \param counts                will be populated such that counts[i][j] will
 *                              be the number of DEs that occurred at
 *                              position j in region/sequence i.
 * \param VERBOSE               if true, output additional status messages to
 *                              stderr about the progress of the computation
 */
void
countDEs(const string &mappedReads_fn, const string &targetRegions_fn,
         const string &experimentType, const string &mapperType,
         vector< vector<size_t> > &counts, const bool VERBOSE = false) {
  const size_t BUFFER_SIZE = 200;
  vector<GenomicRegion> targets;
  ReadBEDFile(targetRegions_fn, targets);

  if (VERBOSE) cerr << "LOADING INTERVAL TREES FROM "
                    << targetRegions_fn << "... ";
  typedef unordered_map<string, IntervalTree<DERoi*, size_t> > TreeMap;
  TreeMap trees;
  loadIntervalTrees_DERoi(targetRegions_fn, trees);
  if (VERBOSE) cerr << "DONE" << endl;

  std::ifstream in(mappedReads_fn.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + mappedReads_fn);

  // for status message
  std::ifstream::pos_type start_of_data = in.tellg();
  in.seekg(0, std::ios::end);
  std::ifstream::pos_type end_of_data = in.tellg();
  in.seekg(start_of_data);
  size_t current_done = 101;

  vector<MappedRead> buffer(BUFFER_SIZE);
  do {
    fill_buffer_mapped_reads(in, mapperType, buffer);
    vector<GenomicRegion> deLocs;
    addDiagEvents(buffer, experimentType, mapperType, deLocs);
    for (size_t i = 0; i < deLocs.size(); ++i) {
      if (trees.count(deLocs[i].get_chrom()) == 0) continue;
      vector<DERoi*> hits;
      TreeMap::iterator t = trees.find(deLocs[i].get_chrom());
      t->second.intersectingInterval(deLocs[i].get_start(),
                                     deLocs[i].get_end(), hits);
      for (size_t j = 0; j < hits.size(); ++j) {
        const size_t relPos = deLocs[i].get_start() - hits[j]->r.get_start();
        if (relPos > hits[j]->counts.size()) {
          stringstream ss;
          ss << "failed counting diagnostic events. DE at location "
             << deLocs[i] << " was found to intersect " << hits[j]->r
             << ", but gave relative position (" << relPos
             << ") beyond range of region";
          throw SMITHLABException(ss.str());
        }
        hits[j]->counts[relPos] += 1;
      }
    }
    size_t percent_done = static_cast<size_t>(in.tellg()) * 100 / end_of_data;
    if ((percent_done != current_done) && (VERBOSE)) {
      std::cerr << "\r" << "PROCESSING " << mappedReads_fn << "... "
                << percent_done << "% completed..." << std::flush;
      current_done = percent_done;
    }
  } while (buffer.size() == BUFFER_SIZE);   // stop when we can't fill
                                            // the buffer anymore (i.e. EOF)
  if (VERBOSE)
    std::cerr << "\r" << "PROCESSING " << mappedReads_fn
              << "... DONE                     " << endl;

  // squash the trees and sort the resultant vector so things have the same
  // order as the input.
  vector<DERoi*> regionsWithCounts;
  for (TreeMap::iterator it = trees.begin(); it != trees.end(); ++it)
    it->second.squash(regionsWithCounts);
  sort(regionsWithCounts.begin(), regionsWithCounts.end(),
      DERoiRankComparator());

  // copy the results and delete the memory allocated by the TreeMap
  for (size_t i = 0; i < regionsWithCounts.size(); ++i) {
    counts.push_back(regionsWithCounts[i]->counts);
  }
  deleteIntervalTrees_DERoi(trees);
}

/**
 * \summary Main entry point for the program; parses command line, extracts
 *          diagnostic events and writes output.
 * \param argc  count of arguments passed to the program
 * \param argv  array of arguments as strings
 */
int
main(int argc, const char **argv) {
  try {
    const string about("This is a program that is part of the Zagros package. "
                       "It's function is to compute the counts of diagnostic "
                       "events within a set of target regions (that the user "
                       "specifies) by extracting them from mapping results.");
    string outputFn = "";
    string regionsFn = "";
    string mapper = "rmap";
    string tech = "iCLIP";
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse(strip_path(argv[0]), about, "<mapped reads>");
    opt_parse.add_opt("output", 'o', "Write output to this file (STDOUT if "
                      "omitted).", OptionParser::OPTIONAL, outputFn);
    opt_parse.add_opt("regions", 'r', "the genomic regions, in BED format, "
                      "corresponding to the input sequences for Zagros.",
                      OptionParser::REQUIRED, regionsFn);
    opt_parse.add_opt("mapper", 'm', "the mapper used to map the reads "
                      "(Default: " + mapper + ")", OptionParser::OPTIONAL,
                      mapper);
    opt_parse.add_opt("tech", 't', "the technology type used in the "
                      "experiment (default " + tech + ")",
                      OptionParser::OPTIONAL, tech);
    opt_parse.add_opt("verbose", 'v', "print more run info (default: " +\
                      toString(VERBOSE), OptionParser::OPTIONAL, VERBOSE);
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

    vector< vector <size_t> > counts;
    countDEs(mappedReadFn, regionsFn, tech, mapper, counts, VERBOSE);
    writeDEs(counts, outputFn, VERBOSE);
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR: " << e.what();
    cerr << endl;
    return EXIT_FAILURE;
  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what();
    cerr << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

