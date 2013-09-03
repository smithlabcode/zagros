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

#ifndef __INPUT_H_
#define __INPUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "ExtendedGenomicRegion.hpp"

struct DE {
  std::string type;
  size_t position;
  std::string base;
  std::string read;
};

struct region_less {
  bool operator()(const GenomicRegion a,
      const GenomicRegion b) const {
    return (a) < (b);
  }
};

class IO {
public:
  /*** Constructors, Destructors and object initialization ***/
  static const size_t flanking_regions_size = 20;

//  static const size_t regions_size = 120;

  static size_t
  adjust_start_pos(const size_t orig_start, const std::string &chrom_name);

  static size_t
  adjust_region_size(const size_t orig_start, const std::string &chrom_name,
      const size_t orig_size);

  static void
  extract_regions_chrom_fasta(const std::string &chrom_name,
      const std::string &filename, const std::vector<GenomicRegion> &regions,
      std::vector<std::string> &sequences, std::vector<std::string> &names);

  static void
  read_dir(const std::string& dirname, std::vector<std::string> &filenames);

  static void
  extract_regions_fasta(const std::string &dirname,
      const std::vector<GenomicRegion> &regions_in,
      std::vector<std::string> &sequences, std::vector<std::string> &names);

  static void
  make_extended_regions_chrom(std::vector<ExtendedGenomicRegion> &regions,
      std::vector<GenomicRegion> &diagnostic_events, const size_t seq_no,
      const std::string experiment, const size_t max_de,
      const size_t min_cluster_size);

  static void
  make_regions_chrom(std::vector<ExtendedGenomicRegion> &regions,
      std::vector<GenomicRegion> &diagnostic_events, const size_t seq_no,
      const std::string experiment, const size_t max_de,
      const size_t min_cluster_size);

  static void
  make_inputs(std::vector<ExtendedGenomicRegion> &mapped_reads,
      std::vector<std::string> &seqs, std::vector<GenomicRegion> &regions,
      std::vector<GenomicRegion> &diagnostic_events,
      const std::string experiment, const size_t max_de,
      const size_t min_cluster_size);

  static void
  make_inputs(std::vector<ExtendedGenomicRegion> &mapped_reads,
      std::vector<GenomicRegion> &regions,
      std::vector<GenomicRegion> &diagnostic_events,
      const std::string experiment, const size_t max_de,
      const size_t min_cluster_size);

  static void
  make_sequence_names(const std::vector<std::string> &names,
      std::vector<std::string> &seqs, const std::vector<GenomicRegion> &regions,
      std::tr1::unordered_map<std::string, size_t> &names_table);

  static void
  fix_sequences(std::vector<ExtendedGenomicRegion> &regions);

  static std::string
  makeAlignment(const std::vector<std::string> &sequences,
      const std::vector<std::vector<double> > &indicators);

  static void
  collapse_extended_regions(std::vector<ExtendedGenomicRegion> &regions,
      std::vector<ExtendedGenomicRegion> &tmp_regions,
      std::vector<GenomicRegion> &diagnostic_events,
      std::vector<std::pair<size_t, std::string> > &boundaries,
      const size_t seq_no, const std::string experiment, const size_t max_de,
      const size_t min_cluster_size);

  static void
  collapse_regions(std::vector<ExtendedGenomicRegion> &regions,
      std::vector<ExtendedGenomicRegion> &tmp_regions,
      std::vector<GenomicRegion> &diagnostic_events,
      std::vector<std::pair<size_t, std::string> > &boundaries,
      const size_t seq_no, const std::string experiment, const size_t max_de,
      const size_t min_cluster_size);

  static void
  load_diagnostic_events(std::vector<GenomicRegion> &de_regions,
      std::tr1::unordered_map<std::string, size_t> &names,
      std::vector<GenomicRegion> &regions,
      std::vector<std::vector<size_t> > &D);

  static void
  parse_diagnostic_events(const std::string &de_string, std::vector<DE> &des);

  static void
  apply_de(ExtendedGenomicRegion &region, DE &de);

  static void
  apply_mutation(ExtendedGenomicRegion &region, DE &de);

  static void
  apply_deletion(ExtendedGenomicRegion &region, DE &de);

  static void
  apply_insertion(ExtendedGenomicRegion &region, DE &de);

  static void
  mergeSequences(ExtendedGenomicRegion &holder, ExtendedGenomicRegion &region);

  static void add_diagnostic_events(
      std::vector<GenomicRegion> &diagnostic_events,
      std::vector<ExtendedGenomicRegion> &containing_reads,
      const std::string name, const std::string experiment,
      const size_t max_de) {
    std::vector<GenomicRegion> tmp_de;
    for (size_t i = 0; i < containing_reads.size(); ++i) {
      add_diagnostic_events(
          tmp_de, containing_reads[i], name, experiment, max_de);
    }
    if (tmp_de.size() > 0) {
      std::tr1::unordered_map<size_t, size_t> number_de;
      for (size_t i = 0; i < tmp_de.size(); i += 1)
        number_de[tmp_de[i].get_start()] += 1;
      for (size_t i = 0; i < tmp_de.size(); i += 1)
        //  if (number_de[tmp_de[i].get_start()] > 100)
        diagnostic_events.push_back(tmp_de[i]);
    }

    /*
     if (tmp_de.size() > 0) {
     std::tr1::unordered_map<size_t,size_t> number_de;
     for ( size_t i = 0; i < tmp_de.size(); i+=1)
     number_de[tmp_de[i].get_start()] += 1;
     size_t max_containing_de = 0;
     size_t max_de_start = 0;
     for ( size_t i = 0; i < tmp_de.size(); i+=1)
     if (number_de[tmp_de[i].get_start()] > max_containing_de) {
     max_containing_de = number_de[tmp_de[i].get_start()];
     max_de_start = tmp_de[i].get_start();
     }
     for ( size_t i = 0; i < tmp_de.size(); i+=1)
     if (tmp_de[i].get_start() == max_de_start)
     diagnostic_events.push_back(tmp_de[i]);
     }
     */
  }

  static void add_diagnostic_events(
      std::vector<GenomicRegion> &diagnostic_events,
      ExtendedGenomicRegion &region, const std::string name,
      const std::string experiment, const size_t max_de) {
    if (diagnostic_events.size() != 0) {
      if (diagnostic_events[diagnostic_events.size() - 1].get_name() != name
          || diagnostic_events[diagnostic_events.size() - 1].get_score()
              < max_de) {
        if (experiment == "iCLIP")
          add_diagnostic_events_iCLIP(diagnostic_events, region, name);
        else if (experiment == "hCLIP")
          add_diagnostic_events_hCLIP(diagnostic_events, region, name);
        else if (experiment == "pCLIP")
          add_diagnostic_events_pCLIP(diagnostic_events, region, name);
      }
    } else {
      if (experiment == "iCLIP")
        add_diagnostic_events_iCLIP(diagnostic_events, region, name);
      else if (experiment == "hCLIP")
        add_diagnostic_events_hCLIP(diagnostic_events, region, name);
      else if (experiment == "pCLIP")
        add_diagnostic_events_pCLIP(diagnostic_events, region, name);
    }
  }

  static void
  add_diagnostic_events_iCLIP(std::vector<GenomicRegion> &diagnostic_events,
      ExtendedGenomicRegion &region, const std::string name);

  static void
  add_diagnostic_events_hCLIP(std::vector<GenomicRegion> &diagnostic_events,
      ExtendedGenomicRegion &region, const std::string name);

  static void
  add_diagnostic_events_pCLIP(std::vector<GenomicRegion> &diagnostic_events,
      ExtendedGenomicRegion &region, const std::string name);

  static void
  find_peak_de_regions(std::vector<GenomicRegion> &de_regions);

  static void
  find_peak_de_regions_chrom(const std::vector<GenomicRegion> &de_regions_chrom,
      std::vector<GenomicRegion> &peak_de_regions);

  static void
  fillTables(const std::vector<std::string> &sequences,
      std::vector<std::vector<double> > &fullStrVectorTable);

  static void
  trimTables(std::vector<std::string> &sequences,
      std::vector<std::vector<double> > &fullStrVectorTable);

  static void
  fillTables(const std::vector<std::string> &sequences,
      std::vector<std::vector<double> > &fullStrVectorTable,
      std::vector<std::vector<std::vector<double> > > &fullStrMatrixTable,
      std::vector<double> &wnSecStr);

  static void
  fillTables(const std::string &fileName,
      const std::vector<std::string> &sequences,
      std::vector<std::vector<double> > &fullStrVectorTable,
      std::vector<std::vector<std::vector<double> > > &fullStrMatrixTable,
      std::vector<double> &wnSecStr);

  static void
  trimTables(std::vector<std::string> &sequences, const std::string &file_name,
      std::vector<std::vector<double> > &fullStrVectorTable,
      std::vector<std::vector<std::vector<double> > > &fullStrMatrixTable);

  static void
  saveTables(const std::vector<std::string> &sequences,
      const std::string &fileName,
      const std::vector<std::vector<double> > &fullStrVectorTable,
      const std::vector<std::vector<std::vector<double> > > &fullStrMatrixTable,
      const std::vector<double> &wnSecStr);

  static void
  save_input_files(const std::vector<std::string> &seqs,
      const std::vector<GenomicRegion> &regions,
      const std::vector<GenomicRegion> &de_regions,
      const std::string base_file);

  static void
  expand_regions(std::vector<GenomicRegion> &regions);

  static void
  unexpand_regions(std::vector<GenomicRegion> &regions);

  static bool
  str_file_checks_out(const GenomicRegion &region,
      const std::string &base_file_name);

  static bool
  str_file_checks_out(const std::string &sequence,
      const std::string &base_file_name);

  static void
  show_percentage(const size_t n, std::string &percentage);

  static void
  read_rmap_output(std::string filename,
      std::vector<ExtendedGenomicRegion> &regions);

  static void
  read_bowtie_output(std::string filename,
      std::vector<ExtendedGenomicRegion> &regions);

  static void
  read_novoalign_output(std::string filename,
      std::vector<ExtendedGenomicRegion> &regions);

  static void
  read_piranha_output(std::string filename,
      std::vector<ExtendedGenomicRegion> &regions);

  static void
  read_piranha_output(std::string filename,
      std::vector<GenomicRegion> &regions);

  static bool
  is_header_line(const std::string& line);

  static bool
  is_track_line(const char *line);

  static void
  convertBowtieExtra(std::string &extra);

  static void
  filter_scores(const float lower_bound, const float upper_bound,
      std::vector<GenomicRegion> &regions);

  static void
  sift_single_chrom(std::vector<GenomicRegion> &other_regions,
      std::vector<GenomicRegion> &regions,
      std::vector<GenomicRegion> &good_regions);

  static void
  sift(std::vector<GenomicRegion> &other_regions,
      std::vector<GenomicRegion> &regions);

private:
  IO();
  /*** Private member constants ***/
};

#endif
