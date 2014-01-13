/**
  \file simulate.cpp
  \brief TODO

  \authors Emad Bahrami Samani, Philip J. Uren

  \section copyright Copyright Details
  Copyright (C) 2012
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


// TODO clean up includes
#include <cstdio>
#include <cstdlib>
#include <sys/time.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "RNG.hpp"
#include "sim_utils.hpp"
#include "OptionParser.hpp"

using std::stringstream;
using std::ofstream;
using std::ostream_iterator;
using std::numeric_limits;
using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::copy;
using std::max;
using std::min;

/*****************************************************************************
 * Constants used by the program
 *****************************************************************************/
namespace SIMRNA {
  const static size_t RANDOM_POSITION = -1;
  const static size_t RNA_ALPHABET_SIZE = 4;
  const static double DEFAULT_BACKGROUND[RNA_ALPHABET_SIZE] = \
                                                      {0.3,  0.2,  0.2,  0.3};

  // types of structure
  const static string ssRNA = "ssRNA";
  const static string dsRNA = "dsRNA";
  const static string unstructured = "unstructured";

  // constants used in creating the structure for the
  // double- and single-stranded motif placement
  const static size_t DEFAULT_OFFSET_FROM_LOOP_START = 3;
  const static size_t DEFAULT_LOOP_LENGTH = 15;
  const static size_t DEFAULT_STEM_LENGTH = 10;

  // DE constants
  const static int MIN_DE_OFFSET = -5;
  const static int MAX_DE_OFFSET = 5;

  inline size_t
  base2int(char c) {
    switch(c) {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    case 'U' : return 3;
    case 'a' : return 0;
    case 'c' : return 1;
    case 'g' : return 2;
    case 't' : return 3;
    case 'u' : return 3;
    default  : return 4;
    }
  }
}

/*****************************************************************************
 * Simple helper functions
 *****************************************************************************/
/***
 * \summary Check s and t for equality, ignoring case. This won't work for
 *          complicated cases, but simple ASCII strings should be fine.
 * \return  true if the string are equal, ignoring case, false otherwise.
 */
bool
isEqualIgnoreCase(string s, string t) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  std::transform(t.begin(), t.end(), t.begin(), ::tolower);
  return s == t;
}

/***
 * \summary Convert a size_t to a string
 */
string
convertSizet(const size_t number) {
  std::stringstream ss;
  ss << number;
  return ss.str();
}

/***
 * \summary Convert an int to a string
 */
string
intToString(const int n) {
  std::stringstream ss;
  ss << n;
  return ss.str();
}

/***
 * \summary check a sorted set of genomic regions to see if any are overlapping
 * \return  true if no regions overlap, false otherwise
 * \throw SMITHLABException if the regions are not sorted
 */
bool
noOverlappingRegions(vector<GenomicRegion> rs) {
 if (!check_sorted(rs))
   throw SMITHLABException("Regions must be in sorted order");
 for (size_t i = 1; i < rs.size(); ++i)
   if (rs[i].overlaps(rs[i-1])) return false;
 return true;
}

static bool
check_consistent(const vector<GenomicRegion> &regions,
    const vector<string> &names, const vector<string> &sequences) {
  if (regions.size() != names.size())
    return false;

  for (size_t i = 0; i < regions.size(); ++i)
    if (regions[i].get_name() != names[i]
        || sequences[i].length() != regions[i].get_width())
      return false;
  return true;
}

/*****************************************************************************
 * Functions for generating random numbers and variables
 *****************************************************************************/

/****
 * @summary: generate a random size_t in a given range (inclusive)
 * @throw:   SMITHLABException if lower is greater than or equal to upper
 */
size_t
randomSize_t(const size_t lower, const size_t upper) {
  if (lower >= upper) {
    stringstream ss;
    ss << "Failed to generate random size_t: "
       << "lower boundary (" << lower << ") was greater than upper ("
       << upper << ").";
      throw SMITHLABException(ss.str());
  }
  return lower + (rand() % (upper - lower + 1));
}


/****
 * @summary: generate a random double in a given range (inclusive)
 * @throw:   SMITHLABException if lower is greater than or equal to upper
 */
double
randomDouble(const double lower, const double upper) {
  if (lower >= upper) {
    stringstream ss;
    ss << "Failed to generate random double: "
       << "lower boundary was greater than upper.";
    throw SMITHLABException(ss.str());
  }
  double frac = (double) rand() / RAND_MAX;
  return lower + (frac * (upper - lower));
}

/***
 * \brief   generate a random int in a given range (inclusive)
 * \throw   SMITHLABException if lower is greater than or equal to upper
 */
int
randomInt(const int lower, const int upper) {
  if (lower >= upper) {
    stringstream ss;
    ss << "Failed to generate random double: "
       << "lower boundary was greater than upper.";
    throw SMITHLABException(ss.str());
  }
  return lower + (rand() % (upper - lower + 1));
}

static size_t
generate_geometric_random_variable(const Runif &rng,
    const double p) {
  return std::floor(std::log(rng.runif(0.0, 1.0)) / std::log(1.0 - p));
}

static vector<double>
generate_dirichlet_random_variables(gsl_rng *gr,
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

static void
generate_targeted_position_weight_matrix(gsl_rng *gr, size_t n,
    size_t m, double target_ic, vector<vector<double> > &pwm) {

  const double pwm_0[SIMRNA::RNA_ALPHABET_SIZE] = { 0.3,  0.2,  0.2,  0.3};
  double ic = 0;
  double sum = 0;
  vector<double> new_col(m, 0.0);
  double ave = 0;
  while (fabs(ave - target_ic) > 0.005) {
    sum = 0;
    pwm.clear();
    while (pwm.size() < n) {
      ic = 0;
      new_col = generate_dirichlet_random_variables(gr, m);
      for (size_t j = 0; j < new_col.size(); ++j)
        if (new_col[j] > 0)
          ic += new_col[j] * fabs(log(new_col[j]/pwm_0[j])/log(2.0));
      sum += ic;
      pwm.push_back(new_col);
    }
    ave = sum/pwm.size();
  }
}

/*****************************************************************************
 * Functions for generating various types of RNA sequences
 *****************************************************************************/

/****
 * @summary: generate a sequence from a nucleotide distribution.
 * @param dist: the distribution to use.
 * @param length: the length of the sequence to build
 * @throw SMITHLABException: if the dist. vector has the wrong dimensions
 */
string
generateSequenceFromNucleotideDistribution(const vector<double> &dist,
                                           const size_t length) {
  if (dist.size() != SIMRNA::RNA_ALPHABET_SIZE) {
    stringstream ss;
    ss << "Failed to generate sequence from nucleotide distribution, "
       << "distribution vector was malformed: found "
       << dist.size() << " entries; expected " << SIMRNA::RNA_ALPHABET_SIZE;
    throw SMITHLABException(ss.str());
  }

  string res = "";
  for (size_t i = 0; i < length; ++i) {
    double r = randomDouble(0, 1);
    if (r < dist[SIMRNA::base2int('A')]) res += 'A';
    else if (r < dist[SIMRNA::base2int('A')] +\
                 dist[SIMRNA::base2int('C')]) res += 'C';
    else if (r < dist[SIMRNA::base2int('A')] +\
                 dist[SIMRNA::base2int('C')] +\
                 dist[SIMRNA::base2int('G')]) res += 'G';
    else res += 'T';
  }

  return res;
}

/****
 * \summary   generate a sequence from a position weight matrix
 * \param pwm the PWM to make the sequence from
 * \throw SMITHLABException if the PWM is malformed.
 */
string
generateSequenceFromPWM(const vector<vector<double> > pwm) {
  string res = "";
  for (size_t j = 0; j < pwm.size(); ++j) {
    if (pwm[j].size() != SIMRNA::RNA_ALPHABET_SIZE) {
      stringstream ss;
      ss << "Failed to generate sequence from PWM, PWM was malformed: found "
         << pwm[j].size() << " entries for position " << j
         << "; expected " << SIMRNA::RNA_ALPHABET_SIZE;
      throw SMITHLABException(ss.str());
    }
    res += generateSequenceFromNucleotideDistribution(pwm[j], 1);
  }
  return res;
}

/*****************************************************************************
 * Functions for placing different types of motif occurrences into sequences
 *****************************************************************************/

/****
 * \brief:     given a PWM, place an occurrence within a sequence with no
 *             particular structure. Note that the motif occurrence will
 *             replace whatever sequence was previously at that place
 * \param seq: the sequence to embed the occurrence in
 * \param pwm: the PWM to generate the occurrence from; should be indexed by
 *             position first, so PWM[i].size() == alphabet size for all i
 * \param indicators    \todo
 * \param pos: location to place the motif occurrence; if equal to
 *             SIMRNA::RANDOM_POSITION, the location is chosen randomly.
 * \throw SMITHLABException if the position is invalid (less than 0, or
 *        greater than seq.size() - motifLength)
 */
size_t
placeUnstructuredMotifOccurrence(string &seq,
                                 const vector<vector<double> > &pwm,
                                 vector<size_t> &indicators,
                                 const size_t pos = SIMRNA::RANDOM_POSITION) {
  // decide where to put the motif
  size_t motifLocation = pos;
  if (motifLocation == SIMRNA::RANDOM_POSITION) {
    if (seq.size() < pwm.size()) {
      stringstream ss;
      ss << "Failed to place unstructured motif occurrence: sequence is "
         << seq.size() << "bp long, but motif is " << pwm.size() << "bp long.";
      throw SMITHLABException(ss.str());
    }
    motifLocation = randomSize_t(0, seq.size() - pwm.size());
  }
  if (motifLocation > seq.size() - pwm.size()) {
    stringstream ss;
    ss << "Failed to place unstructured motif occurrence: invalid "
       << "placement position (" << motifLocation << ") for sequence of "
       << "length " << seq.size();
    throw SMITHLABException(ss.str());
  }

  // place it..
  indicators.push_back(motifLocation);
  string occurrence = generateSequenceFromPWM(pwm);
  seq.replace(seq.begin() + motifLocation,
              seq.begin() + motifLocation + occurrence.size(),
              occurrence);
  return motifLocation;
}


/****
 * \summary given a PWM, place an occurrence within a sequence such that it
 *          falls in a region of particular secondary structure.
 *
 *          The structure is formed by:
 *              if pwm.size() + offset < loopLength
 *              (motif does not extend into second dsRNA sequence)
 *                  take sequence seq[pos-offset-stemLength, pos-offset] and
 *                  place its reverse complement at seq[pos-offset+loopLength]
 *              else (motif extends into second dsRNA sequence, but not first)
 *                  take sequence seq[pos-offset+loopLength] and place reverse
 *                  complement at seq[pos-offset-stemLength, pos-offset]
 *
 *          The occurrence will be 100% in ssRNA if offset > 0 and
 *          pwm.size() + offset < loopLength
 *
 *          X% of the motif will be in dsRNA (from start of motif) if
 *          offset < 0, where x = 100 * min(-offset, pwm.size()) / pwm.size()
 *
 *          X% of the motif will be in dsRNA (from end of motif) if
 *          offset > 0, where
 *          x = 100 * max(0, offset + pwm.size() - loopLength) / pwm.size()
 *
 *
 *          Note that the motif occurrence will replace whatever sequence was
 *          previously at that place. Part of the sequence before/after the
 *          motif will also be modified to force the desired structure to form.
 *
 * \param seq:          the sequence to embed the occurrence in
 * \param pwm:          the PWM to generate the occurrence from; should be
 *                      indexed by position first, so pwm[i].size() == alphabet
 *                      size for all i.
 * \param indicators    \todo
 * \param pos:          location to place the motif occurrence; if equal to
 *                      SIMRNA::RANDOM_POSITION, the location is chosen randomly.
 * \param offset        the distance from the start of the ssRNA loop to the
 *                      start of the motif occurrence. Must satisfy:
 *                      (offset >= -stemLength) and
 *                      (offset <= loopLength + stemLength - pwm.size())
 * \param stemLength    the length of the dsRNA segments intended to form the
 *                      stem
 * \param loopLength    the length of the ssRNA loop at the end of the stem.
 *                      This must be larger than the motif width
 *
 * \throw TODO
 */
size_t
placeStructuredMotifOccurrence(string &seq,
                 const vector<vector<double> > &pwm,
                 vector<size_t> &indicators,
                 const size_t pos = SIMRNA::RANDOM_POSITION,
                 const int offset = SIMRNA::DEFAULT_OFFSET_FROM_LOOP_START,
                 const size_t stemLength = SIMRNA::DEFAULT_STEM_LENGTH,
                 const size_t loopLength = SIMRNA::DEFAULT_LOOP_LENGTH) {
  // check constraints on parameters
  if (loopLength < pwm.size()) {
    stringstream ss;
    ss << "Failed to place structured motif occurrence. Loop length ("
       << loopLength << ") is greater than motif length (" << pwm.size() << ")";
        throw SMITHLABException(ss.str());
  }
  if ((offset < (-1 * static_cast<int>(stemLength))) ||
      (offset > static_cast<int>(loopLength + stemLength - pwm.size()))) {
    stringstream ss;
    ss << "Failed to place structured motif occurrence. Offset ("
       << offset << ") places motif outside of structured region. "
       << "Stem length is " << stemLength << ". Loop length is "
       << loopLength << "." << " Motif length is " << pwm.size() << ". ";
    throw SMITHLABException(ss.str());
  }

  // we need a certain amount of space in the sequence for all this
  // -- complain if the sequence is too short
  if (seq.size() < loopLength + 2 * stemLength) {
    stringstream ss;
    ss << "Failed to place structured motif occurrence. Sequence is too "
       << "short. Need at least " << (loopLength + 2 * stemLength)
       << " nucleotides, but have only " << seq.size() << ".";
    throw SMITHLABException(ss.str());
  }

  // decide where to put the motif
  size_t motifLocation = pos;
  if (motifLocation == SIMRNA::RANDOM_POSITION) {
    // we already checked that the sequence was big enough, so this is safe...
    motifLocation = randomSize_t(stemLength + offset,
                            seq.size() - (stemLength + loopLength - offset));
    indicators.push_back(motifLocation);
  }
  // .. but a user-defined position might not be
  if ((motifLocation < stemLength + offset) ||
      (motifLocation > seq.size() - (stemLength + loopLength - offset))) {
    stringstream ss;
    ss << "Failed to place ssRNA motif occurrence - the chosen start index ("
       << motifLocation << ") is too close to the ";
    if (motifLocation < stemLength + offset) ss << "start ";
    else ss << "end ";
    ss << "of the sequence. There isn't enough space to form the "
       << "secondary structure";
    throw SMITHLABException(ss.str());
  }

  // finally, we can place everything ...
  string occurrence = generateSequenceFromPWM(pwm);
  seq.replace(seq.begin() + motifLocation,
              seq.begin() + motifLocation + occurrence.size(),
              occurrence);
  if (pwm.size() + offset < loopLength) {
    // if the motif doesn't enter the second dsRNA segment, we can safely
    // overwrite it without killing the motif...
    string rcStem = revcomp(seq.substr(motifLocation - offset - stemLength,
                                       stemLength));
    seq.replace(seq.begin() + motifLocation - offset + loopLength,
                seq.begin() + motifLocation - offset + loopLength + stemLength,
                rcStem);
  } else {
    // but if it does, then since pwm.size() < loopLength, we know that it
    // doesn't enter the first dsRNA segment, so we can safely overwrite that
    // instead.
    string rcStem = revcomp(seq.substr(motifLocation - offset + loopLength,
                                       stemLength));
    seq.replace(seq.begin() + motifLocation - offset - stemLength,
                seq.begin() + motifLocation - offset,
                rcStem);
  }
  return motifLocation;
}

/***
 * \summary given a PWM, place an occurrence within a sequence such that it
 *          falls fully within a region of double-stranded RNA
 *
 * \param seq:          the sequence to embed the occurrence in
 * \param pwm:          the PWM to generate the occurrence from; should be
 *                      indexed by position first, so pwm[i].size() == alphabet
 *                      size for all i.
 * \param pos:          location to place the motif occurrence; if equal to
 *                      SIMRNA::RANDOM_POSITION, the location is chosen randomly.
 * \param stemLength    TODO
 * \param loopLength    TODO
 * \return the (rel.) location within the sequence where the motif was placed.
 * \throw TODO
 */
size_t
placeFullDSMotifOccurrence(string &seq,
                 const vector<vector<double> > &pwm,
                 vector<size_t> &indicators,
                 const size_t pos = SIMRNA::RANDOM_POSITION,
                 const size_t stemLength = SIMRNA::DEFAULT_STEM_LENGTH,
                 const size_t loopLength = SIMRNA::DEFAULT_LOOP_LENGTH) {
  // todo PJU: I don't think the below condition is checked anywhere else; the
  // assert is kind of nasty; change this to throw an exception.
  assert(stemLength > pwm.size());

  double coinFlip = (double) rand() / RAND_MAX;
  int oset;
  if (coinFlip > 0.5) {
    // pick an offset that puts the motif in the first dsRNA region
    // must be <= -pwm.size(), but > -stemLength
    oset = randomInt(-stemLength + 1, -pwm.size());
  } else {
    // pick an offset that puts the motif in the second dsRNA region
    // must be > loopLength, but < loopLength + stemLength - pwm.size()
    oset = randomInt(loopLength, loopLength + stemLength - pwm.size());
  }
  return placeStructuredMotifOccurrence(seq, pwm, indicators, pos,
                                        oset, stemLength, loopLength);
}

/****
 * \summary given a PWM, place an occurrence within a sequence such that it
 *          falls fully within a region of single-stranded RNA
 *
 * \param seq:          the sequence to embed the occurrence in
 * \param pwm:          the PWM to generate the occurrence from; should be
 *                      indexed by position first, so pwm[i].size() == alphabet
 *                      size for all i.
 * \param pos:          location to place the motif occurrence; if equal to
 *                      SIMRNA::RANDOM_POSITION, the location is chosen randomly.
 * \param stemLength    TODO
 * \param loopLength    TODO
 * \return the (rel.) location within the sequence where the motif was placed.
 * \throw TODO
 */
size_t
placeFullSSMotifOccurrence(string &seq,
                 const vector<vector<double> > &pwm,
                 vector<size_t> &indicators,
                 const size_t pos = SIMRNA::RANDOM_POSITION,
                 const size_t stemLength = SIMRNA::DEFAULT_STEM_LENGTH,
                 const size_t loopLength = SIMRNA::DEFAULT_LOOP_LENGTH) {
  // pick an offset that puts the motif in the loop. Must be >= 0,
  // but < loopLength - pwm.size()
  double oset = randomSize_t(0, loopLength - pwm.size());
  return placeStructuredMotifOccurrence(seq, pwm, indicators, pos, oset,
                                        stemLength, loopLength);
}

/***
 * \summary place a motif occurrence into the provided sequence using the
 *          provided pwm to generate the occurrence, having the specified
 *          structure type with probability equal to the specified structure
 *          level
 * \param seq        the sequence to place the motif ocurrence into.
 * \param pwm        the pwm to generate the occurrence from.
 * \param indicators the (relative) motif location will be added to this vector;
 *                   any existing entries in this vector are left untouched.
 * \param strType    the structure type to impose on the sequence
 *                   (with prob. equal to strLevel).
 * \param strLevel   0 <= strLevel <= 1; the probability that the motif will
 *                   have the requested structure, instead of being unstructured.
 * \return the (rel.) location within the sequence where the motif was placed.
 * \throw SMITHLABException if <strLevel> isn't between 0 and 1 (inclusive).
 * \throw SMITHLABException if strType does not specify a valid structure type.
 */
size_t
placeMotif(string &seq, const vector< vector<double> > &pwm,
           vector<size_t> &indicators,
           const string &strType = SIMRNA::unstructured,
           const double &strLevel = 1) {
  size_t motifLocation;
  const bool obeyStr = ((randomDouble(0,1) < strLevel));
  if (isEqualIgnoreCase(strType, SIMRNA::unstructured))
    motifLocation = placeUnstructuredMotifOccurrence(seq, pwm, indicators);
  else if (isEqualIgnoreCase(strType, SIMRNA::dsRNA)) {
    if (obeyStr)
      motifLocation = placeFullDSMotifOccurrence(seq, pwm, indicators);
    else
      motifLocation = placeUnstructuredMotifOccurrence(seq, pwm, indicators);
  } else if (isEqualIgnoreCase(strType, SIMRNA::ssRNA)) {
    if (obeyStr)
      motifLocation = placeFullSSMotifOccurrence(seq, pwm, indicators);
    else
      motifLocation = placeUnstructuredMotifOccurrence(seq, pwm, indicators);
  } else {
    throw SMITHLABException(strType + " isn't a valid structure type");
  }
  return motifLocation;
}

/*****************************************************************************
 * Functions for generating diagnostic events
 *****************************************************************************/

/***
 * \summary           place a diagnostic event randomly into the specified
 *                    region such that its distance from a given position is
 *                    geometrically distributed.
 * \param region      the genomic region that the event should occur within.
 * \param diagEvents  the vector of DE locations to place the event into. Any
 *                    existing items in this vector will be left untouched.
 * \param geoP        the parameter P for the geometric distribution to
 *                    determine the distance from the target location;
 *                    0 <= p <= 1
 * \param targetLoc   the (relative to region start) location within the region
 *                    around which DE occurrences will be geometrically
 *                    distributed.
 * \param rng         random number generator to use.
 * \throw SMITHLABException if geoP > 1 or geoP < 0.
 * \throw SMITHLABException if targetLoc is outside of the specified region.
 */
void
placeDEGeom (const GenomicRegion &region, vector<GenomicRegion> &diagEvents,
             const double geoP, const size_t targetLoc,
             const Runif &rng = Runif()) {
  if ((geoP < 0) || (geoP > 1)) {
    stringstream ss;
    ss << "failed to place DE; geometric distribution parameter (" << geoP
       << ") is outside acceptable range (0 <= p <= 1)";
    throw SMITHLABException(ss.str());
  }
  if (region.get_start() + targetLoc > region.get_end()) {
    stringstream ss;
    ss << "failed to place DE; target location (" << targetLoc << ") is "
       << "outside of the specified region " << region;
    throw SMITHLABException(ss.str());
  }

  size_t deLoc;
  do {
    deLoc = generate_geometric_random_variable(rng, geoP);
    if (rng.runif(0, 1) > 0.5) deLoc = region.get_start() + targetLoc + deLoc;
    else deLoc = region.get_start() + targetLoc - deLoc;
  } while ((deLoc < region.get_start()) || (deLoc > region.get_end()));
  diagEvents.push_back(GenomicRegion(region.get_chrom(), deLoc, deLoc + 1,
                                     region.get_name(), 0, region.get_strand()));
}

/***
 * \summary           place a diagnostic event randomly into the specified
 *                    region such that its location follows a random uniform
 *                    distribution.
 * \param region      the genomic region that the event should occur within.
 * \param diagEvents  the vector of DE locations to place the event into. Any
 *                    existing items in this vector will be left untouched.
 * \param rng         random number generator to use.
 */
void
placeDERunif (const GenomicRegion &region, vector<GenomicRegion> &diagEvents,
              const Runif &rng = Runif()) {
  size_t deLoc;
  // need size_t(0) to disambiguate the call
  deLoc = rng.runif(size_t(0), region.get_width());
  diagEvents.push_back(GenomicRegion(region.get_chrom(), deLoc, deLoc + 1,
                                     region.get_name(), 0, region.get_strand()));
}

/***
 * \summary           randomly place a set of diagnostic event into the
 *                    specified region such that their distance from a given
 *                    position follows a geometric distribution
 * \param numDEs      the number of diagnostic events to place
 * \param region      the genomic region that the event should occur within.
 * \param diagEvents  the vector of DE locations to place the event into. Any
 *                    existing items in this vector will be left untouched.
 * \param geoP        the parameter P for the geometric distribution to
 *                    determine the distance from the target location;
 *                    0 <= p <= 1
 * \param targetLoc   the (relative to region start) location within the region
 *                    around which DE occurrences will be geometrically
 *                    distributed.
 * \param rng         random number generator to use.
 */
void
placeDEsGeom (const size_t numDEs, const GenomicRegion &region,
              vector<GenomicRegion> &diagEvents, const double geoP,
              const size_t targetLoc, const Runif &rng = Runif()) {
  for (size_t i = 0; i < numDEs; ++i)
    placeDEGeom(region, diagEvents, geoP, targetLoc, rng);
}

/***
 * \summary           randomly place a set of diagnostic event into the
 *                    specified region such that their locations follow a
 *                    random uniform distribution
 * \param numDEs      the number of diagnostic events to place
 * \param region      the genomic region that the event should occur within.
 * \param diagEvents  the vector of DE locations to place the event into. Any
 *                    existing items in this vector will be left untouched.
 * \param rng         random number generator to use.
 */
void
placeDEsRunif (const size_t numDEs, const GenomicRegion &region,
               vector<GenomicRegion> &diagEvents, const Runif &rng = Runif()) {
  for (size_t i = 0; i < numDEs; ++i)
    placeDERunif(region, diagEvents, rng);
}


/*****************************************************************************
 * IO Functions
 *****************************************************************************/

/****
 * @summary: write the sequences (strings) in <sequences> in fasta format to
 *           <filename> (or STDOUT if <filename> is empty. The sequences are
 *           assumed to correspond to the genomic regions in <regions>,
 *           and sequence[i] will be named using region[i].
 * @throw:   SMITHLABException if (1) the number of regions and sequences
 *           doesn't match or (2) the length of a sequence doesn't match
 *           the length of its corresponding region, or (3) could not open
 *           stream to write to.
 */
void
writeSequencesFasta(const vector<string> &sequences,
                    const vector<GenomicRegion> &regions,
                    const string &filename) {
  // check numbers match up..
  if (sequences.size() != regions.size()) {
    stringstream ss;
    ss << "Failed writing sequences -- number of regions didn't "
       << "match number of sequences";
    throw SMITHLABException(ss.str());
  }

  // make an output stream, use standard out if the output filename is empty.
  ostream *out = (filename.empty()) ? &std::cout :
    new ofstream(filename.c_str());
  if (!out->good()) {
    stringstream ss;
    ss << "Failed to open stream to write to ";
    if (filename.empty()) ss << "STDOUT";
    else ss << filename;
    throw SMITHLABException(ss.str());
  }

  for (size_t i = 0; i < sequences.size(); i++) {
    (*out) << '>' << regions[i].get_name() << '\n' << sequences[i] << '\n';
  }

  if (out != &std::cout) delete out;
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
    for (size_t b = 0; b < SIMRNA::RNA_ALPHABET_SIZE - 1; b++)
      ss << (int) (M[j][b] * N) << "\t";
    ss << (int) (M[j][SIMRNA::RNA_ALPHABET_SIZE - 1] * N) << endl;
  }
  ss << "XX" << endl;
  ss << "AT\tGEO_P=" << p << endl;
  ss << "AT\tGEO_delta=" << delta << endl;
  ss << "XX" << endl;

  for (size_t n = 0; n < N; n++) {
    ss << "BS\t" << sequences[n].substr(indicators[n], M.size()) << "; "
        << regions[n].get_chrom() << ":" << regions[n].get_start() << "-"
        << regions[n].get_end() << "; " << indicators[n] << "; " << M.size()
        << ";  ;" << ((regions[n].get_strand() == '+') ? "p; " : "n;") << endl;
  }
  ss << "XX" << endl;
  ss << "//" << endl;

  return ss.str();
}







/***
 * \summary Count the number of DEs that fall into each of the provided
 *          regions (bins).
 * \param events    must be sorted
 * \param bins      must be non-overlapping and sorted.
 * \param binCounts \todo
 * \throws SMITHLABException \todo
 */
void
countDEsIntoBuckets(const vector<GenomicRegion> &events,
               const vector<GenomicRegion> &bins, vector<size_t> &binCounts) {
  if (!check_sorted(events))
    throw SMITHLABException("DE event locations must be in sorted order");
  if (!check_sorted(bins))
    throw SMITHLABException("Regions must be in sorted order");
  if (!noOverlappingRegions(bins))
    throw SMITHLABException("Regions must be non-overlapping");

  assert(events.size() >= 1);
  assert(bins.size() >= 1);
  binCounts.resize(bins.size(), 0);

  size_t bin_idx = 0, de_idx = 0;
  while (bin_idx < bins.size()) {
    // while we haven't run out of des, and the de chrom is less than the bin,
    // or it's the same but the end of the de is still less than the start of
    // the bin...
    while ((de_idx < events.size()) &&
           ((events[de_idx].get_chrom() < bins[bin_idx].get_chrom()) ||
            ((events[de_idx].get_chrom() == bins[bin_idx].get_chrom()) &&
                (events[de_idx].get_end() < bins[bin_idx].get_start())))) {
      de_idx += 1;
    }
    // if we ran out of des, or the de chrom is greater than the bin, or it's
    // the same but the start of the de is greater than the end of the bin...
    if ((de_idx == events.size()) ||
        (events[de_idx].get_chrom() > bins[bin_idx].get_chrom()) ||
        ((events[de_idx].get_chrom() == bins[bin_idx].get_chrom()) &&
            (events[de_idx].get_start() > bins[bin_idx].get_end())))
      binCounts[bin_idx] = 0;
    else {
      // the de falls inside the bin; count them up while they
      // continue to do so
      while ((de_idx < events.size()) &&
             (bins[bin_idx].overlaps(events[de_idx]))) {
        binCounts[bin_idx] += 1;
        de_idx += 1;
      }
    }
    bin_idx += 1;
  }
}


/****
 * @summary: TODO
 */
int
main(int argc, const char **argv) {
  try {
    string output_basefn;
    size_t deNoiseLevel = 0;
    size_t strLevel = 1;
    string strType = SIMRNA::unstructured;
    double occLevel = 1;
    int site_length = 6;
    double target_ic = 1.7;
    string deFile = "";
    bool VERBOSE = false;
    double geoP = randomDouble(0, 1);
    int geoOffset = randomInt(SIMRNA::MIN_DE_OFFSET, SIMRNA::MAX_DE_OFFSET);

    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse(strip_path(argv[0]),
                                  "program for generating simulated RNAs"
                                  "<fasta-chrom-files>");
    opt_parse.add_opt("output", 'o', "Base name for output files.",
                      OptionParser::REQUIRED, output_basefn);
    opt_parse.add_opt("str_level", 's', "fraction of sequences that "
                      "have the specified structure (0 <= s <= 1). Default: 1",
                      OptionParser::OPTIONAL, strLevel);
    opt_parse.add_opt("occ_level", 'c', "fraction of sequences that have an "
                      "occurrence of the motif. 0 <= c <= 1 (Default: 1)",
                       OptionParser::OPTIONAL, occLevel);
    opt_parse.add_opt("str_type", 't', "type of target site structure (" +\
                      SIMRNA::unstructured + ", " + SIMRNA::ssRNA + ", " +\
                      SIMRNA::dsRNA + "). Default: " + SIMRNA::unstructured,
                      OptionParser::OPTIONAL, strType);
    opt_parse.add_opt("length", 'w', "length of occurrences/PWM (4 <= w "
                      " <= 10) ", OptionParser::OPTIONAL, site_length);
    opt_parse.add_opt("ic", 'i', "target information content of the PWM",
                      OptionParser::OPTIONAL, target_ic);
    opt_parse.add_opt("de_locations", 'd', "bed file with locations of DEs "
                      "from a real experiment; used to simulated number of "
                      "DEs in each sequence", OptionParser::OPTIONAL, deFile);
    opt_parse.add_opt("de_noise", 'n', "fraction of DEs that are noise (i.e. "
                      "not associated with a motif occurrence)",
                      OptionParser::OPTIONAL, deNoiseLevel);
    opt_parse.add_opt("de_p", 'p', "the parameter for the geometric "
                      "distribution governing distance from the motif placement "
                      "to DEs. Default: randomly selected between 0 and 1.",
                      OptionParser::OPTIONAL, geoP);
    opt_parse.add_opt("de_offset", 'f', "offset from the start of motif "
                      "occurrences for the mode of DEs. Default: random "
                      "uniform between " + intToString(SIMRNA::MIN_DE_OFFSET) +\
                      " and " + intToString(SIMRNA::MAX_DE_OFFSET),
                      OptionParser::OPTIONAL, geoOffset);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    if ((site_length < 4 || site_length > 10)) {
      cerr << "target length should satisfy: (3 < target length <= 10)"
           << endl;
      return EXIT_SUCCESS;
    }
    const string input_bed_file(leftover_args.front());
    const string input_fasta_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    // this doesn't make sense, so don't let the user do it
    if ((strLevel != 1) && (isEqualIgnoreCase(strType, SIMRNA::unstructured))) {
      throw SMITHLABException("structure level parameter cannot be used when "
                              "simulating unstructured motifs");
    }

    // check up front that a valid structure was selected, so we don't waste
    // time if it's junk.
    if ((!isEqualIgnoreCase(strType, SIMRNA::unstructured)) &&
        (!isEqualIgnoreCase(strType, SIMRNA::dsRNA)) &&
        (!isEqualIgnoreCase(strType, SIMRNA::ssRNA)))
      throw SMITHLABException(strType + " is not a valid structure type");

    // the first step is to read the regions and the sequences the user has
    // given us, making sure that that they match up, are sorted and
    // non-overlapping.
    vector<GenomicRegion> regions;
    vector<string> seqs, names;
    if (VERBOSE) cerr << "READING REGIONS FILE... ";
    ReadBEDFile(input_bed_file, regions);
    if (VERBOSE) cerr << "DONE (READ " << regions.size() << " REGIONS)"
                      << endl << "READING SEQUENCES FILE... ";
    read_fasta_file(input_fasta_file, names, seqs);
    if (VERBOSE) cerr << "DONE (READ " << seqs.size() << " SEQUENCES)"
                      << endl << "CHECKING CONSISTENCY... ";
    if (!check_consistent(regions, names, seqs))
      throw SMITHLABException("inconsistent region names in FASTA/BED files");
    if (VERBOSE) cerr << "DONE (FILES MATCH)" << endl;

    // seed the random number generator using the time, and make a uniform
    // random number generator object, and a gsl random number generator
    if (VERBOSE) cerr << "MAKING RANDOM NUMBER GENERATOR... ";
    struct timeval time;
    gettimeofday(&time,NULL);
    size_t seed = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    srand(seed);
    const Runif rng(seed);
    gsl_rng *gr = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gr, static_cast<unsigned long int >(seed));
    if (VERBOSE) cerr << "DONE" << endl;

    // read the distribution of the number of diagnostic events per sequence
    vector<size_t> numDEsPerSeq;
    if (!deFile.empty()) {
      if (VERBOSE) cerr << "READING DIAGNOSTIC EVENTS... ";
      vector<GenomicRegion> deLocations;
      ReadBEDFile(deFile, deLocations);
      countDEsIntoBuckets(deLocations, regions, numDEsPerSeq);
      if (VERBOSE) cerr << "DONE (FOUND " << deLocations.size()
                        << " DIAG. EVENTS)" << endl;
    } else {
      numDEsPerSeq.resize(seqs.size(), 0);
    }

    // generate the PWM for the motif we're going to place
    if (VERBOSE) cerr << "GENERATING PWM... ";
    vector<vector<double> > pwm;
    generate_targeted_position_weight_matrix(gr, site_length,
                                             SIMRNA::RNA_ALPHABET_SIZE,
                                             target_ic, pwm);
    if (VERBOSE) {
      cerr << "DONE. GENERATED MATRIX: " << endl;
      for (size_t j = 0; j < pwm.size(); j++) {
        for (size_t k = 0; k < pwm[j].size(); k++)
          cerr << pwm[j][k] << ", ";
        cerr << endl;
      }
    }

    // figure out what the output filenames will be
    string outfile_seq = output_basefn + ".fa";
    string outfile_mat = output_basefn + ".mat";
    string outfile_de = output_basefn + "_de.bed";

    // place the motif occurrences and diagnostic events into the sequences
    vector<size_t> indicators;
    vector<GenomicRegion> diagnosticEvents;
    if (VERBOSE) cerr << "SIMULATING MOTIF OCCURRENCES... ";
    for (size_t i = 0; i < seqs.size(); ++i) {
      // first, decide whether we will even place a motif in this sequence
      bool place = (randomDouble(0,1) > occLevel);
      size_t numDEsToPlace = numDEsPerSeq[i], numRealDEsPlaced = 0;

      if (place) {
        // place the motif with the desired structure
        const size_t motifLoc = placeMotif(seqs[i], pwm, indicators,
                                           strType, strLevel);

        // place the real DEs
        assert(numDEsToPlace >= 0);
        // todo check that deNoiseLevel is >= 0 and <= 1
        numRealDEsPlaced = ceil(numDEsToPlace * deNoiseLevel);
        placeDEsGeom(numRealDEsPlaced, regions[i], diagnosticEvents, geoP,
                     motifLoc);
      }

      // places the noise DEs for this sequence
      placeDEsRunif(numDEsToPlace-numRealDEsPlaced,
                    regions[i], diagnosticEvents);
    }
    if (VERBOSE) cerr << "DONE." << endl;

    // time to write output...
    if (VERBOSE) cerr << "WRITING OUTPUT... ";
    // sequences
    writeSequencesFasta(seqs, regions, outfile_seq);
    // matrix
    ofstream out_mat(outfile_mat.c_str());
    out_mat << print_model(pwm, geoP, geoOffset, regions, seqs, indicators);
    // diagnostic events
    ofstream de_out(outfile_de.c_str());
    copy(diagnosticEvents.begin(), diagnosticEvents.end(),
         std::ostream_iterator<GenomicRegion>(de_out, "\n"));
    if (VERBOSE) cerr << "DONE" << endl;;

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
