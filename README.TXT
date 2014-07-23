
                    ______                         
                   |___  /                         
                      / / __ _  __ _ _ __ ___  ___ 
                     / / / _` |/ _` | '__/ _ \/ __|
                    / /_| (_| | (_| | | | (_) \__ \
                   /_____\__,_|\__, |_|  \___/|___/
                                __/ |              
                               |___/               
                   ********************************
                   *            V1.1.0            *
                   ********************************


*********************************
Copyright and License Information
*******************************************************************************
Copyright (C) 2013
University of Southern California,
Emad Bahrami-Samani, Philip J. Uren, Andrew D. Smith
  
Authors: Emad Bahrami-Samani, Philip J. Uren, Andrew D. Smith
  
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


***********************
Building and Installing 
*******************************************************************************
This software package has been designed to operate in a UNIX-like environment.
It has been tested on MacOS X Snow Leopard and Linux. 

Step 0
------
  This software package requires a functioning installation of the GNU 
  Scientific Library (GSL). If you don't already have this installed, you 
  will need to download and install it from http://www.gnu.org/software/gsl/

Step 1
------
  To build the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution  
  
  > make all 
  
Step 2
------
  To install the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution
  
  > make install
  
  This will place the binaries in the bin directory under the package root.
  They can be used directly from there without any additional steps. You can
  add that directory to your PATH environment variable to avoid having to 
  specify their full paths, or you can copy the binaries to another directory
  of your choice in your PATH 
  
  
***********
Basic Usage
*******************************************************************************

Zagros has four modes of operations.


******************
1.) Sequence only: 

In this mode, only the sequence information is used for motif discovery. The 
input can either be a set of sequences in fasta format or genomic regions in 
bed format. This set of regions/sequences corresponds to the locations of 
significant enrichment for reads in the experiment.

In case the input is the set of sequences, you can simply run:

    > ./zagros input.fa

If the input consists of a set of genomic regions, the set of the target genome
sequences must also be provided to Zagros to extract the sequnces:

    > ./zagros -c path/to/chrom_directory input.bed

The chromosome directory can be downloaded from UCSC genome browser website.
The lastest versions can be found here:
http://hgdownload.soe.ucsc.edu/downloads.html


**************************
2.) Sequence and Structure

In this mode, in addition to the target sequences the secondary structure information
is used as well. In this case, the secondary structure data must be first obtained 
and saved using the "thermo" program:

    > ./thermo -o input.str input.fa

    or

    > ./thermo -c path/to/chrom_directory -o input.str input.bed

After this step, by providing both the target and secondary structure file to Zagros 
the motif discovery is performed based on both.

    > ./zagros -t input.str input.fa

    or

    > ./zagros -c path/to/chrom_directory -t input.str input.bed


**********************************
3.) Sequence and Diagnostic events

In this mode, in addition to the target sequences the information about cross-link
modification events is used as well. In this case, the diagnostic events information 
must be first obtained and saved using the "extractDEs" program.

The input to extractDEs program is the set of mapped reads. The user must specify
what technology is used for obtaining the reads (hCLIP, pCLIP or iCLIP), what mapper
is used for mapping the reads, and the genomic regions of significant regions that
is used for zagros as input. ExtractDEs then produces the set of diagnostic events 
corrsponding to the regions of interest. Zagros can interprete the mapped reads from 
three mappers: bowtie (native output format), novoalign (native output format) and 
RMAP (bed format).

    > ./extractDEs -m novoalign -t iCLIP -o input.des -r input.bed mapped_reads.novo

Then run zagros program by inputing the diagnostic events file as one of the options.

    > ./zagros -d input.des input.fa

    or

    > ./zagros -c path/to/chrom_directory -d input.des input.bed


*********************************************
4.) Sequence, Structure and Diagnostic events
    > ./thermo -o input.str input.fa
    > ./extractDEs -m novoalign -t iCLIP -o input.des -r mapped_reads.nov
    > ./zagros -t input.str -d input.des input.fa


********************
Command Line Options
*******************************************************************************

Usage: zagros [OPTIONS] <target_regions/sequences>

Options:
  -o, -output             output file name (default: stdout) 
  -w, -width              width of motifs to find (4 <= w <= 12; default: 6) 
  -n, -number             number of motifs to output (default: 1) 
  -c, -chrom              directory with chrom files (FASTA format) 
  -t, -structure          structure information file 
  -d, -diagnostic_events  diagnostic events information file 
  -i, -diagEventsThresh   down-sample diagnostic events to this many per 
                          sequence (-1 for no down-sampling; default: -1) 
  -s, -starting-points    number of starting points to try for EM search. Higher 
                          values will be slower, but more likely to find the 
                          global maximum. 
  -v, -verbose            print more run info 

Help options:
  -?, -help               print this help message 
      -about              print about message 


**************************************
4. Input and output formats of Zagros
*******************************************************************************

Zagros takes three possible inputs for mapped reads:
  1.) RMAP      (BED format)              [REQUIRED if using diagnostic evnts]
  2.) Novoalign (Novoalign native format) [REQUIRED if using diagnostic evnts]      
  3.) Bowtie    (Bowtie native format)    [REQUIRED if using diagnostic evnts]


************************
Contacts and bug reports
*******************************************************************************
Emad Bahrami-Samani 
bahramis@usc.edu

Philip J. Uren
uren@usc.edu

Andrew D. Smith
andrewds@usc.edu

If you found a bug in Piranha, we'd like to know about it. Before contacting us
though, please check the following list:

  1.) Are you using the latest version? The bug you found may already have 
      been fixed.
  2.) Check that your input is in the correct format and you have selected
      the correct options.
  3.) Please reduce your input to the smallest possible size that still 
      produces the bug; we will need your input data to reproduce the 
      problem.


