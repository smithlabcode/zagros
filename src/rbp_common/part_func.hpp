/**
 \file part_func.hpp
 \brief TODO

 \todo what other license info and origin acknowledgments need to be given
 for this file?

 \section copyright Copyright Details
 Copyright (C) 2012
 University of Southern California,
 Emad Bahrami Samani, Philip J. Uren, Andrew D. Smith

 \authors Emad Bahrami Samani, Philip J. Uren

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
 The classes in this file don't follow proper const-correctness yet

 \section revisions Revision Details

 ****/

#ifndef __PARTFUNC_H_
#define __PARTFUNC_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

/** \todo see if these are used anywhere and replace them with STL max/min **/
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))

/** \todo replace these with constants **/
#define NST 0
#define DEF -50
#define NSM 0
#define GASCONST 1.98717
#define K0  273.15
#define INF 1000000
#define FORBIDDEN 9999
#define BONUS 10000
#define NBPAIRS 7
#define TURN 3
#define MAXLOOP 30
#define ENCODE(c) encode_char(c)
#define NBASES 8
#define SCALE 10
#define MAXALPHA 20
/** \brief TODO **/
#define temperature 37.0
/** \brief temperature of param measurements */
#define Tmeasure 37+K0
/** \brief parameter for logarithmic loop */
#define lxc37 107.856

/** \todo yuck yuck yuck **/
#define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):\
	SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))

#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif
#ifdef WITH_DMALLOC
#define space(S) calloc(1,(S))
#endif
#undef xrealloc

/**
 * \brief A class for calculating base pair probabilities for RNAs using
 *        J.S. McCaskill's partition function algorithm.
 */
class MC {
public:
  /*** Constructors, destructors and object initialisation ***/
  MC();
  ~MC();
  double getTemperature();
  void init_pf_fold(int length, float energy);
  void free_pf_arrays(void);
  void getProbVector(std::vector<double> &vec);
  void getProbMatrix(std::vector<std::vector<double> > &mat, int size);
  float getMinimumFreeEnergy(char *sequence, char *structure);
  float pf_fold(char *sequence, char *structure);
  double expLoopEnergy(int u1, int u2, int type, int type2, short a1, short a2,
      short a3, short a4);

private:
  /*** Private Member Functions ***/
  void *space(unsigned int size);
  void nrerror(const char message[]);
  void update_pf_params(int length, float energy);
  char bppm_symbol(float *x);
  void printMatrix(std::string fname, int n, int seql, int pos);
  double expHairpinEnergy(int u, int type, short si1, short sj1,
      const char *string);
  void sprintf_bppm(int length, char *structure);
  void scale_pf_params(unsigned int length, float energy);
  void get_arrays(unsigned int length);
  void make_ptypes(const short *S, const char *structure);
  void pf_linear(char *sequence, char *structure);
  void pf_create_bppm(char *sequence, char *structure);
  int encode_char(char c);
  void make_pair_matrix(void);
  void* xrealloc(void *p, unsigned size);
  std::string convertInt(int number);

  /************************ Private Member Variables *************************/
  double *pr; /* base pairing prob. matrix */
  int *iindx; /* pr[i,j] -> pr[iindx[i]-j] */
  double expMLclosing, expMLintern[NBPAIRS + 1], *expMLbase;
  double expTermAU;
  double expdangle5[NBPAIRS + 1][5], expdangle3[NBPAIRS + 1][5];
  double lxc, exptetra[40], expTriloop[40];
  double expstack[NBPAIRS + 1][NBPAIRS + 1];
  double expmismatchI[NBPAIRS + 1][5][5], expmismatchH[NBPAIRS + 1][5][5],
      expmismatchM[NBPAIRS + 1][5][5];
  double expint11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  double expint21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  double expint22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  double *exphairpin;
  double expbulge[MAXLOOP + 1];
  double expinternal[MAXLOOP + 1];
  double expninio[5][MAXLOOP + 1];

  double *q, *qb, *qm, *qm1, *qqm, *qqm1, *qq, *qq1;
  double *prml, *prm_l, *prm_l1, *q1k, *qln;
  double *scale;
  char *ptype; /* precomputed array of pair types */
  int *qtype; /* precomputed array of pair types */
  int *mqtype; /* precomputed array of pair types */
  int *jindx;
  char *sequence;

  short alias[MAXALPHA + 1];
  int pair[MAXALPHA + 1][MAXALPHA + 1];

  short *S, *S1;

  int mismatchM37[NBPAIRS + 1][5][5];

  /** \brief scaling factor to avoid floating point overflows */
  double pf_scale;

  /** \bug this is not set anywhere, but it is used ... **/
  int Triloop_E37[40];

  /************************ Private Static Constants *************************/
  /** \brief constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*k + ML_BASE*u **/
  const static int ML_BASE37 = 0;
  /** \brief constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*k + ML_BASE*u **/
  const static int ML_closing37 = 340;
  /** \brief constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*k + ML_BASE*u **/
  const static int ML_intern37 = 40;

  /** \brief max Ninio-correction for asymmetric internal loops with branches
   *         n1 and n2 \f$ ninio_energy =
   *         min{max_ninio, |n1-n2|*F_ninio[min{4.0,n1,n2}] } \f$
   *  \todo  fix this equation..
   */
  const static int MAX_NINIO = 300;
  /* only F[2] used */

  /** \brief stabilizing contribution due to special hairpins of size 4
   (tetraloops) */
  const static int TETRA_ENTH37 = -400;

  /** \brief penalty for AU (or GU) terminating helix. Mismatches already
   contain these */
  const static int TerminalAU = 50;

  /** \brief penalty for forming a bi-molecular duplex */
  const static int DuplexInit = 410;

};

#endif
