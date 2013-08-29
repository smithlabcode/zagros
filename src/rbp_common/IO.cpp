/*    sortbed: a program for sorting BED format files
 *    Copyright (C) 2008 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include <numeric>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "ExtendedGenomicRegion.hpp"
#include "Model.hpp"
#include "IO.hpp"
#include "RNA_Utils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::accumulate;
using std::ostream;
using std::pair;
using std::make_pair;
using std::sort;
using std::numeric_limits;

using smithlab::alphabet_size;

using std::tr1::unordered_map;

typedef ExtendedGenomicRegion* ExtendedGenomicRegionPointer;

IO::IO() {
}

void IO::make_regions_chrom(vector<ExtendedGenomicRegion> &regions,
    vector<GenomicRegion> &diagnostic_events, const size_t seq_no,
    const string experiment, const size_t max_de,
    const size_t min_cluster_size) {

  vector<pair<size_t, string> > pos_boundaries;
  vector<pair<size_t, string> > neg_boundaries;

  for (size_t i = 0; i < regions.size(); ++i) {
    if (regions[i].pos_strand()) {
      string s = "f" + convertSizet(i);
      pos_boundaries.push_back(make_pair(regions[i].get_start() - 1, s));
      s = "t" + convertSizet(i);
      pos_boundaries.push_back(make_pair(regions[i].get_end() - 1, s));
    } else {
      string s = "f" + convertSizet(i);
      neg_boundaries.push_back(make_pair(regions[i].get_start() - 1, s));
      s = "t" + convertSizet(i);
      neg_boundaries.push_back(make_pair(regions[i].get_end() - 1, s));
    }
  }
  vector<ExtendedGenomicRegion> tmp_regions = regions;
  regions.clear();

  collapse_regions(
      regions, tmp_regions, diagnostic_events, pos_boundaries, seq_no,
      experiment, max_de, min_cluster_size);
  collapse_regions(
      regions, tmp_regions, diagnostic_events, neg_boundaries, seq_no,
      experiment, max_de, min_cluster_size);

}

void IO::make_extended_regions_chrom(vector<ExtendedGenomicRegion> &regions,
    vector<GenomicRegion> &diagnostic_events, const size_t seq_no,
    const string experiment, const size_t max_de,
    const size_t min_cluster_size) {

  vector<pair<size_t, string> > pos_boundaries;
  vector<pair<size_t, string> > neg_boundaries;

  for (size_t i = 0; i < regions.size(); ++i) {
    if (regions[i].pos_strand()) {
      string s = "f" + convertSizet(i);
      pos_boundaries.push_back(make_pair(regions[i].get_start() - 1, s));
      s = "t" + convertSizet(i);
      pos_boundaries.push_back(make_pair(regions[i].get_end() - 1, s));
    } else {
      string s = "f" + convertSizet(i);
      neg_boundaries.push_back(make_pair(regions[i].get_start() - 1, s));
      s = "t" + convertSizet(i);
      neg_boundaries.push_back(make_pair(regions[i].get_end() - 1, s));
    }
  }
  vector<ExtendedGenomicRegion> tmp_regions = regions;
  regions.clear();

  collapse_extended_regions(
      regions, tmp_regions, diagnostic_events, pos_boundaries, seq_no,
      experiment, max_de, min_cluster_size);
  collapse_extended_regions(
      regions, tmp_regions, diagnostic_events, neg_boundaries, seq_no,
      experiment, max_de, min_cluster_size);

}

void IO::sort_regions(vector<GenomicRegion> &regions, string outfile) {
  vector<GenomicRegionPointer> sorter;
  for (vector<GenomicRegion>::iterator i = regions.begin(); i != regions.end();
      ++i)
    sorter.push_back(&(*i));
  sort(sorter.begin(), sorter.end(), region_pointer_less());
  ostream* out =
      (outfile.empty()) ? &std::cout : new std::ofstream(outfile.c_str());
  for (vector<GenomicRegionPointer>::const_iterator i(sorter.begin());
      i != sorter.end(); ++i)
    *out << *(*i) << '\n';
  if (out != &std::cout)
    delete out;
}

size_t IO::adjust_start_pos(const size_t orig_start, const string &chrom_name) {
  static const double LINE_WIDTH = 50.0;
  const size_t name_offset = chrom_name.length() + 2; // For the '>' and the '\n';
  const size_t preceding_newlines = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  return orig_start + preceding_newlines + name_offset;
}

size_t IO::adjust_region_size(const size_t orig_start, const string &chrom_name,
    const size_t orig_size) {
  static const double LINE_WIDTH = 50.0;
  const size_t preceding_newlines_start = static_cast<size_t>(std::floor(
      orig_start / LINE_WIDTH));
  const size_t preceding_newlines_end = static_cast<size_t>(std::floor(
      (orig_start + orig_size) / LINE_WIDTH));
  return (orig_size + (preceding_newlines_end - preceding_newlines_start));
}

void IO::extract_regions_chrom_fasta(const string &chrom_name,
    const string &filename, const vector<GenomicRegion> &regions,
    vector<string> &sequences, vector<string> &names) {

  std::ifstream in(filename.c_str());
  for (vector<GenomicRegion>::const_iterator i(regions.begin());
      i != regions.end(); ++i) {

    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;

    const size_t start_pos = adjust_start_pos(orig_start_pos, chrom_name);
    const size_t region_size = adjust_region_size(
        orig_start_pos, chrom_name, orig_region_size);
    assert(start_pos >= 0);

    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);

    std::remove_if(
        buffer, buffer + region_size,
        std::bind2nd(std::equal_to<char>(), '\n'));
    buffer[orig_region_size] = '\0';

    sequences.push_back(buffer);
    names.push_back(i->get_name());
    std::transform(
        sequences.back().begin(), sequences.back().end(),
        sequences.back().begin(), std::ptr_fun(&toupper));
    if (i->neg_strand())
      revcomp_inplace(sequences.back());
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}

void IO::read_dir(const string& dirname, vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw "could not open directory: " + dirname;

  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  if (errno)
    throw "error reading directory: " + dirname;
  if (filenames.empty())
    throw "no valid files found in: " + dirname;
  closedir(dir);
}

void IO::extract_regions_fasta(const string &dirname,
    const vector<GenomicRegion> &regions_in, vector<string> &sequences,
    vector<string> &names) {

  static const string FASTA_SUFFIX(".fa");
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<GenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);

  std::tr1::unordered_map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;

  for (size_t i = 0; i < regions.size(); ++i) {
    const string chrom_name(regions[i].front().get_chrom());
    const string chrom_file(chrom_name + FASTA_SUFFIX);
    std::tr1::unordered_map<string, size_t>::const_iterator f_idx =
        chrom_regions_map.find(chrom_file);
    if (f_idx == chrom_regions_map.end())
      throw SMITHLABException("chrom not found:\t" + chrom_file);
    extract_regions_chrom_fasta(
        chrom_name, filenames[f_idx->second], regions[i], sequences, names);
  }
}

void IO::collapse_regions(vector<ExtendedGenomicRegion> &regions,
    vector<ExtendedGenomicRegion> &tmp_regions,
    vector<GenomicRegion> &diagnostic_events,
    vector<pair<size_t, string> > &boundaries, const size_t seq_no,
    const string experiment, const size_t max_de,
    const size_t min_cluster_size) {

  static const string FAKE_NAME("X");
  const string chrom(regions.front().get_chrom());

  ExtendedGenomicRegion holder(chrom, 0, 0, FAKE_NAME, 0, '+', "", "");
  GenomicRegion de_holder(chrom, 0, 0, FAKE_NAME, 0, '+');
  sort(boundaries.begin(), boundaries.end());

  vector<ExtendedGenomicRegion> containing_reads;

  size_t count = 0;
  size_t score = 0;
  for (size_t i = 0; i < boundaries.size(); ++i) {
    size_t index = convertString(
        (boundaries[i].second).substr(1, (boundaries[i].second).length() - 1));
    if ((boundaries[i].second).substr(0, 1) == "t") {
      --count;
      if (count == 0) {
        holder.set_end(boundaries[i].first);
        holder.set_score(score);
        if (score > min_cluster_size) {
          regions.push_back(holder);
          add_diagnostic_events(
              diagnostic_events, containing_reads, holder.get_name(),
              experiment, max_de);
        }
        score = 0;
      }
    } else {
      if (count == 0) {
        containing_reads.clear();
        holder.set_start(boundaries[i].first);
        holder.set_end(tmp_regions[index].get_end());
        string name = "sequence_" + convertSizet(seq_no + regions.size() + 1);
        holder.set_name(name);
        holder.set_strand(tmp_regions[index].get_strand());
        containing_reads.push_back(tmp_regions[index]);
        ++count;
        ++score;
      } else if (abs(holder.get_start() - boundaries[i].first) > 500) {
        holder.set_end(boundaries[i].first);
        holder.set_score(score);
        if (score > min_cluster_size) {
          regions.push_back(holder);
          add_diagnostic_events(
              diagnostic_events, containing_reads, holder.get_name(),
              experiment, max_de);
        }
        holder.set_start(boundaries[i].first);
        holder.set_end(tmp_regions[index].get_end());
        string name = "sequence_" + convertSizet(seq_no + regions.size() + 1);
        holder.set_name(name);
        holder.set_strand(tmp_regions[index].get_strand());
        containing_reads.push_back(tmp_regions[index]);
        count = 1;
        score = 1;
      } else {
        if (holder.get_start() < tmp_regions[index].get_start()
            && holder.get_end() < tmp_regions[index].get_end()
            && holder.get_end() >= tmp_regions[index].get_start()) {
          holder.set_end(tmp_regions[index].get_end());
        }
        containing_reads.push_back(tmp_regions[index]);
        ++count;
        ++score;
      }
    }
  }
}

void IO::collapse_extended_regions(vector<ExtendedGenomicRegion> &regions,
    vector<ExtendedGenomicRegion> &tmp_regions,
    vector<GenomicRegion> &diagnostic_events,
    vector<pair<size_t, string> > &boundaries, const size_t seq_no,
    const string experiment, const size_t max_de,
    const size_t min_cluster_size) {

  static const string FAKE_NAME("X");
  const string chrom(regions.front().get_chrom());

  ExtendedGenomicRegion holder(chrom, 0, 0, FAKE_NAME, 0, '+', "", "");
  GenomicRegion de_holder(chrom, 0, 0, FAKE_NAME, 0, '+');
  sort(boundaries.begin(), boundaries.end());

  vector<ExtendedGenomicRegion> containing_reads;

  size_t count = 0;
  size_t score = 0;
  for (size_t i = 0; i < boundaries.size(); ++i) {
    size_t index = convertString(
        (boundaries[i].second).substr(1, (boundaries[i].second).length() - 1));
    if ((boundaries[i].second).substr(0, 1) == "t") {
      --count;
      if (count == 0) {
        holder.set_end(boundaries[i].first);
        holder.set_score(score);
        if (score > min_cluster_size) {
          regions.push_back(holder);
          add_diagnostic_events(
              diagnostic_events, containing_reads, holder.get_name(),
              experiment, max_de);
        }
        score = 0;
      }
    } else {
      if (count == 0) {
        containing_reads.clear();
        holder.set_start(boundaries[i].first);
        holder.set_end(tmp_regions[index].get_end());
        string name = "sequence_" + convertSizet(seq_no + regions.size() + 1);
        holder.set_name(name);
        holder.set_strand(tmp_regions[index].get_strand());
        holder.set_sequence(tmp_regions[index].get_sequence());
        containing_reads.push_back(tmp_regions[index]);
        ++count;
        ++score;
      } else if (abs(holder.get_start() - boundaries[i].first) > 500) {
        holder.set_end(boundaries[i].first);
        holder.set_score(score);
        if (score > min_cluster_size) {
          regions.push_back(holder);
          add_diagnostic_events(
              diagnostic_events, containing_reads, holder.get_name(),
              experiment, max_de);
        }
        holder.set_start(boundaries[i].first);
        holder.set_end(tmp_regions[index].get_end());
        string name = "sequence_" + convertSizet(seq_no + regions.size() + 1);
        holder.set_name(name);
        holder.set_strand(tmp_regions[index].get_strand());
        holder.set_sequence(tmp_regions[index].get_sequence());
        containing_reads.push_back(tmp_regions[index]);
        count = 1;
        score = 1;
      } else {
        if (holder.get_start() < tmp_regions[index].get_start()
            && holder.get_end() < tmp_regions[index].get_end()
            && holder.get_end() >= tmp_regions[index].get_start()) {
          mergeSequences(holder, tmp_regions[index]);
          holder.set_end(tmp_regions[index].get_end());
        }
        containing_reads.push_back(tmp_regions[index]);
        ++count;
        ++score;
      }
    }
  }
}

void IO::add_diagnostic_events_iCLIP(vector<GenomicRegion> &diagnostic_events,
    ExtendedGenomicRegion &region, const string name) {

  size_t score = 1;
  if (diagnostic_events.size() != 0)
    if (diagnostic_events[diagnostic_events.size() - 1].get_name() == name)
      score = diagnostic_events[diagnostic_events.size() - 1].get_score() + 1;
  if (region.pos_strand()) {
    GenomicRegion gr(
        region.get_chrom(), region.get_start() - 1, region.get_start(), name,
        score, region.get_strand());
    diagnostic_events.push_back(gr);
  } else {
    GenomicRegion gr(
        region.get_chrom(), region.get_start() - 2, region.get_start() - 1,
        name, score, region.get_strand());
    diagnostic_events.push_back(gr);
  }
}

void IO::add_diagnostic_events_hCLIP(vector<GenomicRegion> &diagnostic_events,
    ExtendedGenomicRegion &region, const string name) {

  vector<DE> des;
  parse_diagnostic_events(region.get_extra(), des);
  if (des.size() == 1) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events[diagnostic_events.size() - 1].get_name() == name)
        score = diagnostic_events[diagnostic_events.size() - 1].get_score() + 1;
    for (size_t j = 0; j < des.size(); ++j) {
      if (region.pos_strand() && des[j].type == "insertion"
          && des[j].base == "T") {
        GenomicRegion gr(
            region.get_chrom(), region.get_start() + des[j].position - 1,
            region.get_start() + des[j].position, name, score,
            region.get_strand());
        diagnostic_events.push_back(gr);
      } else if (!region.pos_strand() && des[j].type == "insertion"
          && des[j].base == "A") {
        GenomicRegion gr(
            region.get_chrom(), region.get_start() + des[j].position - 2,
            region.get_start() + des[j].position - 1, name, score,
            region.get_strand());
        diagnostic_events.push_back(gr);
      }
    }
  }
}

void IO::add_diagnostic_events_pCLIP(vector<GenomicRegion> &diagnostic_events,
    ExtendedGenomicRegion &region, const string name) {

  vector<DE> des;
  parse_diagnostic_events(region.get_extra(), des);
  if (des.size() == 1) {
    size_t score = 1;
    if (diagnostic_events.size() != 0)
      if (diagnostic_events[diagnostic_events.size() - 1].get_name() == name)
        score = diagnostic_events[diagnostic_events.size() - 1].get_score() + 1;
    for (size_t j = 0; j < des.size(); ++j) {
      if (region.pos_strand() && des[j].type == "mutation" && des[j].base == "T"
          && des[j].read == "C") {
        GenomicRegion gr(
            region.get_chrom(), region.get_start() + des[j].position - 1,
            region.get_start() + des[j].position, name, score,
            region.get_strand());
        diagnostic_events.push_back(gr);
      } else if (!region.pos_strand() && des[j].type == "mutation"
          && des[j].base == "A" && des[j].read == "G") {
        GenomicRegion gr(
            region.get_chrom(), region.get_start() + des[j].position - 2,
            region.get_start() + des[j].position - 1, name, score,
            region.get_strand());
        diagnostic_events.push_back(gr);
      }
    }
  }
}

void IO::make_inputs(vector<ExtendedGenomicRegion> &mapped_reads,
    vector<string> &seqs, vector<GenomicRegion> &regions,
    vector<GenomicRegion> &diagnostic_events, const string experiment,
    const size_t max_de, const size_t min_cluster_size) {

  vector<vector<ExtendedGenomicRegion> > separated_by_chrom;
  separate_extended_chromosomes(mapped_reads, separated_by_chrom);
  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    fix_sequences(separated_by_chrom[i]);
    make_extended_regions_chrom(
        separated_by_chrom[i], diagnostic_events, regions.size(), experiment,
        max_de, min_cluster_size);
    for (size_t j = 0; j < separated_by_chrom[i].size(); ++j) {
      GenomicRegion gr(
          separated_by_chrom[i][j].get_chrom(),
          separated_by_chrom[i][j].get_start(),
          separated_by_chrom[i][j].get_end(),
          separated_by_chrom[i][j].get_name(),
          separated_by_chrom[i][j].get_score(),
          separated_by_chrom[i][j].get_strand());
      regions.push_back(gr);
      seqs.push_back(separated_by_chrom[i][j].get_sequence());
    }
    separated_by_chrom[i].clear();
  }
}

void IO::make_inputs(vector<ExtendedGenomicRegion> &mapped_reads,
    vector<GenomicRegion> &regions, vector<GenomicRegion> &diagnostic_events,
    const string experiment, const size_t max_de,
    const size_t min_cluster_size) {

  vector<vector<ExtendedGenomicRegion> > separated_by_chrom;
  separate_extended_chromosomes(mapped_reads, separated_by_chrom);
  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    make_regions_chrom(
        separated_by_chrom[i], diagnostic_events, regions.size(), experiment,
        max_de, min_cluster_size);
    for (size_t j = 0; j < separated_by_chrom[i].size(); ++j) {
      GenomicRegion gr(
          separated_by_chrom[i][j].get_chrom(),
          separated_by_chrom[i][j].get_start(),
          separated_by_chrom[i][j].get_end(),
          separated_by_chrom[i][j].get_name(),
          separated_by_chrom[i][j].get_score(),
          separated_by_chrom[i][j].get_strand());
      regions.push_back(gr);
      cerr << "\r"
          << int(
              (100 * j * i)
                  / (separated_by_chrom.size() * separated_by_chrom[i].size()))
          << "% completed..." << std::flush;
    }
    separated_by_chrom[i].clear();
  }
  cerr << "\r" << "100% completed..." << endl;
}

void IO::fix_sequences(vector<ExtendedGenomicRegion> &regions) {

  for (size_t i = 0; i < regions.size(); ++i) {
    vector<DE> des;
    parse_diagnostic_events(regions[i].get_extra(), des);
    for (size_t j = 0; j < des.size(); ++j) {
      apply_de(regions[i], des[j]);
    }
  }
}

void IO::mergeSequences(ExtendedGenomicRegion &holder,
    ExtendedGenomicRegion &region) {

  if (region.get_end() > holder.get_end()
      && holder.get_end() >= region.get_start()) {
    if (holder.pos_strand()) {
      holder.set_sequence(
          holder.get_sequence()
              + region.get_sequence().substr(
                  holder.get_end() - region.get_start(),
                  region.get_end() - holder.get_end()));
    } else {
      if (holder.get_end() == region.get_start())
        holder.set_sequence(holder.get_sequence() + region.get_sequence());
      else
        holder.set_sequence(
            region.get_sequence()
                + holder.get_sequence().substr(
                    holder.get_end() - region.get_start(),
                    holder.get_sequence().length()
                        - abs(holder.get_end() - region.get_start())));
    }
  }
}

void IO::apply_de(ExtendedGenomicRegion &region, DE &de) {
  if (de.type == "mutation")
    apply_mutation(region, de);
  else if (de.type == "deletion")
    apply_deletion(region, de);
  else if (de.type == "insertion")
    apply_insertion(region, de);
}

void IO::apply_mutation(ExtendedGenomicRegion &region, DE &de) {
  if (region.pos_strand())
    region.set_sequence(
        region.get_sequence().replace(de.position - 1, 1, de.base));
  else
    region.set_sequence(
        region.get_sequence().replace(
            region.get_sequence().length() - de.position, 1,
            RNAUtils::getComplement(de.base)));
}

void IO::apply_deletion(ExtendedGenomicRegion &region, DE &de) {
  if (region.pos_strand())
    region.set_sequence(
        region.get_sequence().erase(de.position - 1, de.base.length()));
  else
    region.set_sequence(
        region.get_sequence().erase(
            region.get_sequence().length() - (de.position - 1)
                - de.base.length(), de.base.length()));
  region.set_end(region.get_start() + region.get_sequence().length());
}

void IO::apply_insertion(ExtendedGenomicRegion &region, DE &de) {
  if (region.pos_strand())
    region.set_sequence(region.get_sequence().insert(de.position - 1, de.base));
  else
    region.set_sequence(
        region.get_sequence().insert(
            region.get_sequence().length() - (de.position - 1),
            RNAUtils::getReverseComplement(de.base)));
  region.set_end(region.get_start() + region.get_sequence().length());
}

void IO::parse_diagnostic_events(const string &de_string, vector<DE> &des) {
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
      } else {
        index = parts[i].find("+");
        if (index != string::npos) {
          de.type = "deletion";
          de.position = convertString(parts[i].substr(0, index));
          de.base = parts[i].substr(index + 1, parts[i].length() - index - 1);
        } else {
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

void IO::make_sequence_names(const vector<string> &names, vector<string> &seqs,
    const vector<GenomicRegion> &regions,
    unordered_map<string, size_t> &names_table) {

  vector<string> sorted_seqs;
  for (size_t i = 0; i < regions.size(); ++i) {
    for (size_t j = 0; j < names.size(); ++j)
      if (regions[i].get_name() == names[j]) {
        sorted_seqs.push_back(seqs[j]);
        break;
      }
    names_table[regions[i].get_name()] = i;
  }
  seqs.clear();
  seqs = sorted_seqs;
}

void IO::find_peak_de_regions(vector<GenomicRegion> &de_regions) {

  vector<vector<GenomicRegion> > de_regions_chrom;
  separate_chromosomes(de_regions, de_regions_chrom);
  vector<GenomicRegion> peak_de_regions;
  for (size_t i = 0; i < de_regions_chrom.size(); ++i)
    find_peak_de_regions_chrom(de_regions_chrom[i], peak_de_regions);
  de_regions.clear();
  de_regions = peak_de_regions;
}

void IO::find_peak_de_regions_chrom(
    const vector<GenomicRegion> &de_regions_chrom,
    vector<GenomicRegion> &peak_de_regions) {

  if (de_regions_chrom.size() > 0) {
    std::tr1::unordered_map<size_t, size_t> number_de;
    std::tr1::unordered_map<string, size_t> name_de;

    for (size_t i = 0; i < de_regions_chrom.size(); i += 1) {
      number_de[de_regions_chrom[i].get_start()] = 0;
      name_de[de_regions_chrom[i].get_name()] = 0;
    }

    for (size_t i = 0; i < de_regions_chrom.size(); i += 1) {
      number_de[de_regions_chrom[i].get_start()] += 1;
    }

    for (size_t i = 0; i < de_regions_chrom.size(); i += 1)
      if (number_de[de_regions_chrom[i].get_start()] > 0)
        if (number_de[de_regions_chrom[i].get_start()]
            > name_de[de_regions_chrom[i].get_name()])
          name_de[de_regions_chrom[i].get_name()] =
              number_de[de_regions_chrom[i].get_start()];

    for (size_t i = 0; i < de_regions_chrom.size(); i += 1)
      if (number_de[de_regions_chrom[i].get_start()] > 0)
        if (number_de[de_regions_chrom[i].get_start()]
            == name_de[de_regions_chrom[i].get_name()]) {
          peak_de_regions.push_back(de_regions_chrom[i]);
          peak_de_regions.back().set_score(
              number_de[de_regions_chrom[i].get_start()]);
          number_de[de_regions_chrom[i].get_start()] = 0;
          name_de[de_regions_chrom[i].get_name()] = 0;
        }
  }
}

void IO::load_diagnostic_events(vector<GenomicRegion> &de_regions,
    unordered_map<string, size_t> &names, vector<GenomicRegion> &regions,
    vector<vector<size_t> > &D) {
  D.resize(regions.size());
  for (size_t i = 0; i < de_regions.size(); ++i)
    if (regions[names[de_regions[i].get_name()]].pos_strand())
      D[names[de_regions[i].get_name()]].push_back(
          de_regions[i].get_start()
              - regions[names[de_regions[i].get_name()]].get_start()
              + IO::flanking_regions_size - 1);
    else
      D[names[de_regions[i].get_name()]].push_back(
          regions[names[de_regions[i].get_name()]].get_width()
              + 2 * (IO::flanking_regions_size - 1)
              - (de_regions[i].get_start()
                  - regions[names[de_regions[i].get_name()]].get_start()
                  + IO::flanking_regions_size - 1));
}

string IO::makeAlignment(const vector<string> &sequences,
    const vector<vector<double> > &indicators) {
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

  double max_X = -1;
  int max_i = -1;
  stringstream ss;
  for (size_t n = 0; n < N; n++) {
    max_X = -1;
    max_i = -1;
    for (size_t i = 0; i < indicators[n].size(); i++) {
      if (indicators[n][i] > max_X) {
        max_X = indicators[n][i];
        max_i = i;
      }
    }
    ss << max_i + 1 << "\t"
        << sequences[n].substr(
            max_i, sequences[n].length() - indicators[n].size() + 1) << endl;
  }
  return ss.str();
}

void IO::fillTables(const vector<string> &sequences,
    vector<vector<double> > &fullStrVectorTable) {

  size_t N = sequences.size();

  cerr << "Calculating the structure information..." << endl;
  fullStrVectorTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrVectorTable[n].resize(L);
  }
  for (size_t n = 0; n < N; n++) {
    cerr << "\r" << int((100 * n) / N) << "% completed..." << std::flush;
    fullStrVectorTable[n] = RNAUtils::getBasePairProbabilityVector(
        sequences[n]);
  }
  cerr << "\r" << "100% completed..." << endl;
}

void IO::trimTables(vector<string> &sequences,
    vector<vector<double> > &fullStrVectorTable) {

  vector<vector<double> > exp_fullStrVectorTable = fullStrVectorTable;

  size_t N = sequences.size();

  vector<string> trimmed_sequences;
  for (size_t n = 0; n < N; n++)
    trimmed_sequences.push_back(
        sequences[n].substr(
            flanking_regions_size,
            sequences[n].length() - (2 * flanking_regions_size)));
  sequences.clear();
  sequences = trimmed_sequences;

  fullStrVectorTable.clear();
  fullStrVectorTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrVectorTable[n].resize(L);
    for (size_t i = 0; i < L; i++)
      fullStrVectorTable[n][i] = exp_fullStrVectorTable[n][i
          + flanking_regions_size];
  }
}

void IO::fillTables(const vector<string> &sequences,
    vector<vector<double> > &fullStrVectorTable,
    vector<vector<vector<double> > > &fullStrMatrixTable,
    vector<double> &wnSecStr) {

  size_t N = sequences.size();

  cerr << "Calculating the structure information..." << endl;
  fullStrVectorTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrVectorTable[n].resize(L);
  }

  fullStrMatrixTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrMatrixTable[n].resize(L);
    for (size_t i = 0; i < L; i++)
      fullStrMatrixTable[n][i].resize(L);
  }

  vector<vector<vector<double> > > weightedNumberSecStrTable;
  weightedNumberSecStrTable.resize(N);

  vector<bool> large_flag(N, false);
  for (size_t n = 0; n < N; n++) {
    cerr << "\r" << int((100 * n) / N) << "% completed..." << std::flush;
    size_t L = sequences[n].length();
    fullStrVectorTable[n] = RNAUtils::getBasePairProbabilityVector(
        sequences[n]);
    double sum_l = 0;

    for (size_t l = 0; l < sequences[n].length(); l++)
      sum_l += fullStrVectorTable[n][l];
    if (sum_l <= 0 || sum_l > L)
      large_flag[n] = true;
    if (large_flag[n]) {
      fullStrVectorTable[n].clear();
      fullStrVectorTable[n].resize(L, 0.5);
    }

    if (!large_flag[n])
      fullStrMatrixTable[n] = RNAUtils::getBasePairProbabilityMatrix(
          sequences[n]);
    else {
      fullStrMatrixTable[n].clear();
      fullStrMatrixTable[n].resize(L);
      for (size_t i = 0; i < L; i++)
        fullStrMatrixTable[n][i].resize(L, 1 / L);
    }

    if (!large_flag[n])
      weightedNumberSecStrTable[n] =
          RNAUtils::calculateWeightedNumberOfStructures(
              sequences[n], fullStrMatrixTable[n]);
    else {
      weightedNumberSecStrTable[n].clear();
      weightedNumberSecStrTable[n].resize(L);
      for (size_t i = 0; i < L; i++)
        weightedNumberSecStrTable[n][i].resize(L, 1);
    }

    wnSecStr.push_back(
        weightedNumberSecStrTable[n][flanking_regions_size][L
            - flanking_regions_size - 1]);
  }

  cerr << "\r" << "100% completed..." << endl;
}

void IO::fillTables(const string &fileName, const vector<string> &sequences,
    vector<vector<double> > &fullStrVectorTable,
    vector<vector<vector<double> > > &fullStrMatrixTable,
    vector<double> &wnSecStr) {
  cerr << "Structure information available - Loading data...";
  string inputFileName;
  ifstream inf;

  inputFileName = fileName + ".str";
  inf.open(inputFileName.c_str());
  size_t N = sequences.size();

  wnSecStr.clear();
  for (size_t n = 0; n < N; n++) {
    string s;
    getline(inf, s);
    double d;
    stringstream s2d(s);
    s2d >> d;
    wnSecStr.push_back(d);
  }
  inf.close();

  fullStrVectorTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrVectorTable[n].resize(L);
  }

  vector<bool> large_flag(N, false);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrVectorTable[n] = RNAUtils::getBasePairProbabilityVector(
        sequences[n]);
    double sum_l = 0;

    for (size_t l = 0; l < sequences[n].length(); l++)
      sum_l += fullStrVectorTable[n][l];
    if (sum_l <= 0 || sum_l > L)
      large_flag[n] = true;
    if (large_flag[n]) {
      fullStrVectorTable[n].clear();
      fullStrVectorTable[n].resize(L, 0.5);
    }
  }

  fullStrMatrixTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrMatrixTable[n].resize(L);
    for (size_t i = 0; i < L; i++)
      fullStrMatrixTable[n][i].resize(L);
  }
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    if (!large_flag[n])
      fullStrMatrixTable[n] = RNAUtils::getBasePairProbabilityMatrix(
          sequences[n]);
    else {
      fullStrMatrixTable[n].clear();
      fullStrMatrixTable[n].resize(L);
      for (size_t i = 0; i < L; i++)
        fullStrMatrixTable[n][i].resize(L, 1 / L);
    }
  }
  cerr << "done!" << endl;
}

void IO::trimTables(vector<string> &sequences, const string &file_name,
    vector<vector<double> > &fullStrVectorTable,
    vector<vector<vector<double> > > &fullStrMatrixTable) {

  vector<vector<double> > exp_fullStrVectorTable = fullStrVectorTable;
  vector<vector<vector<double> > > exp_fullStrMatrixTable = fullStrMatrixTable;

  size_t N = sequences.size();

  vector<string> trimmed_sequences;
  for (size_t n = 0; n < N; n++)
    trimmed_sequences.push_back(
        sequences[n].substr(
            flanking_regions_size,
            sequences[n].length() - (2 * flanking_regions_size)));
  sequences.clear();
  sequences = trimmed_sequences;

  fullStrVectorTable.clear();
  fullStrVectorTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrVectorTable[n].resize(L);
    for (size_t i = 0; i < L; i++)
      fullStrVectorTable[n][i] = exp_fullStrVectorTable[n][i
          + flanking_regions_size];
  }

  fullStrMatrixTable.clear();
  fullStrMatrixTable.resize(N);
  for (size_t n = 0; n < N; n++) {
    size_t L = sequences[n].length();
    fullStrMatrixTable[n].resize(L);
    for (size_t i = 0; i < L; i++) {
      fullStrMatrixTable[n][i].resize(L);
      for (size_t j = 0; j < L; j++)
        fullStrMatrixTable[n][i][j] = exp_fullStrMatrixTable[n][i
            + flanking_regions_size][j + flanking_regions_size];
    }
  }
}

void IO::saveTables(const vector<string> &sequences, const string &fileName,
    const vector<vector<double> > &fullStrVectorTable,
    const vector<vector<vector<double> > > &fullStrMatrixTable,
    const vector<double> &wnSecStr) {
  string outputFileName;
  ofstream outf;

  size_t N = sequences.size();

  outputFileName = fileName + ".str";
  outf.open(outputFileName.c_str());

  for (size_t n = 0; n < N; n++)
    outf << wnSecStr[n] << endl;
  outf.close();
}

void IO::save_input_files(const vector<string> &seqs,
    const vector<GenomicRegion> &regions,
    const vector<GenomicRegion> &de_regions, const string base_file) {

  string seqs_outfile = base_file + ".fa";
  std::ostream* out1 =
      (!seqs_outfile.empty()) ? new ofstream(seqs_outfile.c_str()) : &cout;
  for (size_t i = 0; i < seqs.size(); ++i) {
    *out1 << ">" << regions[i].get_name() << endl << seqs[i] << endl;
  }
  if (out1 != &cout)
    delete out1;

  string de_outfile = base_file + "_de.bed";
  std::ostream* out3 =
      (!de_outfile.empty()) ? new ofstream(de_outfile.c_str()) : &cout;
  copy(
      de_regions.begin(), de_regions.end(),
      std::ostream_iterator<GenomicRegion>(*out3, "\n"));
  if (out3 != &cout)
    delete out3;
}

void IO::expand_regions(vector<GenomicRegion> &regions) {

  for (size_t i = 0; i < regions.size(); ++i) {
    size_t width = regions[i].get_width();
    if (regions[i].get_width() <= regions_size) {
      regions[i].set_start(
          regions[i].get_start()
              - int((regions_size - width) / 2)
              - flanking_regions_size);
      regions[i].set_end(
          regions[i].get_end()
              + int((regions_size - width) / 2)
              + flanking_regions_size);
    } else {
      regions[i].set_start(
          regions[i].get_start()
              + int((width - regions_size) / 2)
              - flanking_regions_size);
      regions[i].set_end(
          regions[i].get_end()
              - int((width - regions_size) / 2)
              + flanking_regions_size);
    }
  }
//  for (size_t i = 0; i < regions.size(); ++i) {
//    regions[i].set_start(regions[i].get_start() - flanking_regions_size);
//    regions[i].set_end(regions[i].get_end() + flanking_regions_size);
//  }
}

void IO::unexpand_regions(vector<GenomicRegion> &regions) {
  for (size_t i = 0; i < regions.size(); ++i) {
    regions[i].set_start(regions[i].get_start() + flanking_regions_size);
    regions[i].set_end(regions[i].get_end() - flanking_regions_size);
  }
}

bool IO::str_file_checks_out(const GenomicRegion &region,
    const string &base_file_name) {

  string file_name = base_file_name + ".str";
  ifstream in(file_name.c_str());
  if (!in)
    return false;
  else
    return true;
}

bool IO::str_file_checks_out(const string &sequence,
    const string &base_file_name) {

  string file_name = base_file_name + ".str";
  ifstream in(file_name.c_str());
  if (!in)
    return false;
  else
    return true;
}

void IO::show_percentage(const size_t n, string &percentage) {
  percentage = convertSizet(n) + "% [";
  for (size_t i = 0; i < n - 1; ++i)
    percentage = percentage + "=";
  for (size_t i = n; i < 100; ++i)
    percentage = percentage + " ";
  percentage = percentage + "]";
}

void IO::read_rmap_output(string filename,
    vector<ExtendedGenomicRegion> &regions) {
  static const size_t buffer_size = 10000;

  std::ifstream in(filename.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + filename);

  std::ifstream::pos_type start_of_data = in.tellg();
  in.seekg(0, std::ios::end);
  std::ifstream::pos_type end_of_data = in.tellg();
  in.seekg(start_of_data);

  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    if (!is_header_line(buffer) && !is_track_line(buffer)) {
      vector<string> parts(split_extended_whitespace_quoted(buffer));
      char strand = ((parts[5] == "+") ? '+' : '-');
      ExtendedGenomicRegion gr(
          parts[0], convertString(parts[1]), convertString(parts[2]), parts[3],
          convertString(parts[4]), strand, parts[6], "2C>C");
      regions.push_back(gr);
    }

    size_t percent_done = static_cast<size_t>(in.tellg()) * 100 / end_of_data;
    std::cerr << "\r" << percent_done << "% completed..." << std::flush;
    in.peek();
  }
  in.close();
  std::cerr << std::endl;
}

void IO::read_bowtie_output(string filename,
    vector<ExtendedGenomicRegion> &regions) {
  static const size_t buffer_size = 10000;

  std::ifstream in(filename.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + filename);

  std::ifstream::pos_type start_of_data = in.tellg();
  in.seekg(0, std::ios::end);
  std::ifstream::pos_type end_of_data = in.tellg();
  in.seekg(start_of_data);

  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    if (!is_header_line(buffer) && !is_track_line(buffer)) {
      vector<string> parts(split_extended_whitespace_quoted(buffer));
      if (parts.size() == 8) {
        vector<string> name(smithlab::split_whitespace_quoted(parts[0]));
        char strand = ((parts[1] == "+") ? '+' : '-');
        convertBowtieExtra(parts[7]);
        ExtendedGenomicRegion gr(
            parts[2], convertString(parts[3]),
            convertString(parts[3]) + parts[4].length(), name[0],
            convertString(parts[6]), strand, parts[4], parts[7]);
        regions.push_back(gr);
      }
    }
    size_t percent_done = static_cast<size_t>(in.tellg()) * 100 / end_of_data;
    std::cerr << "\r" << percent_done << "% completed..." << std::flush;
    in.peek();
  }
  in.close();
  std::cerr << std::endl;
}

void IO::read_novoalign_output(string filename,
    vector<ExtendedGenomicRegion> &regions) {
  static const size_t buffer_size = 10000;

  std::ifstream in(filename.c_str());
  if (!in)
    throw BEDFileException("cannot open input file " + filename);

  std::ifstream::pos_type start_of_data = in.tellg();
  in.seekg(0, std::ios::end);
  std::ifstream::pos_type end_of_data = in.tellg();
  in.seekg(start_of_data);

  while (!in.eof()) {
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw BEDFileException("Line too long in file: " + filename);
    if (!is_header_line(buffer)) {
      vector<string> parts(split_extended_whitespace_quoted(buffer));
      if (parts.size() == 14)
        if (parts[4] == "U") {
          vector<string> name(smithlab::split_whitespace_quoted(parts[0]));
          char strand = ((parts[9] == "F") ? '+' : '-');
          ExtendedGenomicRegion gr(
              parts[7].substr(1), convertString(parts[8]),
              convertString(parts[8]) + parts[2].length(), name[0],
              convertString(parts[5]), strand, parts[2], parts[13]);
          regions.push_back(gr);
        }
    }
    size_t percent_done = static_cast<size_t>(in.tellg()) * 100 / end_of_data;
    std::cerr << "\r" << percent_done << "% completed..." << std::flush;
    in.peek();
  }
  in.close();
  std::cerr << std::endl;
}

void IO::read_piranha_output(string filename,
    vector<ExtendedGenomicRegion> &regions) {

  vector<GenomicRegion> the_regions;
  ReadBEDFile(filename, the_regions);
  for (size_t i = 0; i < the_regions.size(); i += 1) {
    if (the_regions[i].get_width() <= regions_size) {
      ExtendedGenomicRegion gr(
          the_regions[i].get_chrom(),
          the_regions[i].get_start(),
          the_regions[i].get_end(),
          the_regions[i].get_name(),
          the_regions[i].get_score(),
          the_regions[i].get_strand(), "ACGT", "2C>C");
      regions.push_back(gr);
    }
  }
}

void IO::read_piranha_output(string filename, vector<GenomicRegion> &regions) {

  vector<GenomicRegion> the_regions;
  ReadBEDFile(filename, the_regions);
  for (size_t i = 0; i < the_regions.size(); i += 1) {
      string name = "sequence_" + convertSizet(i+1);
      GenomicRegion gr(
          the_regions[i].get_chrom(),
          the_regions[i].get_start(),
          the_regions[i].get_end(),
          name,
          the_regions[i].get_score(),
          the_regions[i].get_strand());
      regions.push_back(gr);
  }
}

bool IO::is_header_line(const string& line) {
  if (line.substr(0, 1) == "#")
    return true;
  static const char *browser_label = "browser";
  static const size_t browser_label_len = 7;
  for (size_t i = 0; i < browser_label_len; ++i)
    if (line[i] != browser_label[i])
      return false;
  return true;
}

bool IO::is_track_line(const char *line) {
  static const char *track_label = "track";
  static const size_t track_label_len = 5;
  for (size_t i = 0; i < track_label_len; ++i)
    if (line[i] != track_label[i])
      return false;
  return true;
}

void IO::convertBowtieExtra(string &extra) {
  std::replace(extra.begin(), extra.end(), ',', ' ');
  extra.erase(std::remove(extra.begin(), extra.end(), ':'), extra.end());
}

void IO::filter_scores(const float lower_bound, const float upper_bound,
    vector<GenomicRegion> &regions) {
  vector<GenomicRegion> new_regions;
  new_regions.reserve(regions.size() / 2);
  for (size_t i = 0; i < regions.size(); ++i) {
    const double score(regions[i].get_score());
    if (score >= lower_bound && score <= upper_bound)
      new_regions.push_back(regions[i]);
  }
  regions.swap(new_regions);
}

void IO::sift_single_chrom(vector<GenomicRegion> &other_regions,
    vector<GenomicRegion> &regions, vector<GenomicRegion> &good_regions) {

  typedef vector<GenomicRegion>::iterator region_itr;

  for (size_t i = 0; i < regions.size(); ++i) {
    region_itr closest(find_closest(other_regions, regions[i]));
    if (closest->overlaps(regions[i])) {
      good_regions.push_back(regions[i]);
      good_regions[good_regions.size() - 1].set_name(closest->get_name());
    }
  }
}

void IO::sift(vector<GenomicRegion> &other_regions,
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
    const unordered_map<string, size_t>::const_iterator j = chrom_lookup.find(
        regions_by_chrom[i].front().get_chrom());
    if (j != chrom_lookup.end()) {
      sift_single_chrom(
          other_regions_by_chrom[j->second], regions_by_chrom[i], good_regions);
    }

  }
  regions.swap(good_regions);


  for (size_t i = 0; i < other_regions.size(); ++i) {
    for (size_t j = 0; j < regions.size(); ++j)
      if (other_regions[i].get_name() == regions[j].get_name()) {
        other_regions[i].set_strand(regions[j].get_strand());
        break;
      }
  }
}

