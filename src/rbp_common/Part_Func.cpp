/**
 \file part_func.cpp
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

 \section revisions Revision Details

 ****/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <float.h> 
#include <iostream>
#include <fstream>
#include <bitset>
#include <limits>

#include "Part_Func.hpp"
#include "RNA_Utils.hpp"
#include "Int_Loops.hpp"
//#include "expLoopTable.h"

using std::string;
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::endl;

/******************************************************************************
 * Constructors and destructors for the MC class
 */

/**
 * \brief TODO
 * \todo fix magic number -1
 */
MC::MC() {
  pf_scale = -1;
  S = NULL;
  S1 = NULL;
}

/****
 * \brief cleanup dynamically allocated memory this MC object is using
 */
MC::~MC() {
  free_pf_arrays();
}

/****
 * \brief TODO
 */
int MC::encode_char(char c) {
  int code;
  const char *pos;
  pos = strchr(Law_and_Order, c);
  if (pos == NULL)
    code = 0;
  else
    code = (int) (pos - Law_and_Order);
  if (code > 4)
    code--; /* make T and U equivalent */
  return code;
}

/****
 * \brief TODO
 */
void MC::make_pair_matrix(void) {
  int i, j;

  for (i = 0; i < 5; i++)
    alias[i] = (short) i;
  alias[5] = 3; /* X <-> G */
  alias[6] = 2; /* K <-> C */
  alias[7] = 0; /* I <-> default base '@' */
  for (i = 0; i < NBASES; i++) {
    for (j = 0; j < NBASES; j++)
      pair[i][j] = BP_pair[i][j];
  }
  pair[3][4] = pair[4][3] = 0;
  for (i = 0; i < NBASES; i++) {
    for (j = 0; j < NBASES; j++)
      rtype[pair[i][j]] = pair[j][i];
  }
}

/****
 * \brief TODO
 */
void* MC::space(unsigned size) {
  void *pointer;

  if ((pointer = (void *) calloc(1, (size_t) size)) == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"SPACE: requested size: %d\n", size);
      nrerror("SPACE allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
    nrerror("SPACE allocation failure -> no memory");
  }
  return pointer;
}

/****
 * \brief TODO
 */
void* MC::xrealloc(void *p, unsigned size) {
  if (p == 0)
    return space(size);
  p = (void *) realloc(p, size);
  if (p == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"xrealloc: requested size: %d\n", size);
      nrerror("xrealloc allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
    nrerror("xrealloc allocation failure -> no memory");
  }
  return p;
}

/****
 * \brief TODO
 */
void MC::nrerror(const char message[]) /* output message upon error */
{
  fprintf(stderr, "\n%s\n", message);
  exit(EXIT_FAILURE);
}

/**
 * \brief populate the given vector with the probability matrix from this
 *        MC object.
 * \param size the length of the sequence that was folded.
 */
void MC::getProbMatrix(std::vector<std::vector<double> > &mat, size_t seq_length) {

  std::vector<double> vec;
  for (size_t i = 1; i <= seq_length; i++) {
    for (size_t j = i; j <= seq_length; j++) {
      vec.push_back(pr[iindx[i] - j]);
    }
    mat.push_back(vec);
  }
}

/**
 * \brief TODO
 */
float MC::pf_fold(char *sequence, char *structure) {
  double Q;

  double free_energy;
  int n = (int) strlen(sequence);

  /* do the linear pf fold and fill all matrices  */
  pf_linear(sequence, structure);

  Q = q[iindx[1] - n];

  /* ensemble free energy in Kcal/mol              */
  if (Q <= FLT_MIN) {
    fprintf(stderr, "\npf_scale too large\n");
    cout << sequence << endl;
    cout << structure << endl;
//	cout << (-log(Q)-n*log(pf_scale))*(temperature+K0)*GASCONST/1000.0 << endl;
//        exit (1);
  }

  free_energy = (-log(Q) - n * log(pf_scale)) * (temperature + K0) * GASCONST
      / 1000.0;
  /* in case we abort because of floating point errors */
  if (n > 1600)
    fprintf(stderr, "free energy = %8.2f\n", free_energy);

  /* calculate base pairing probability matrix (bppm)  */
  pf_create_bppm(sequence, structure);

  return Q;
}

float MC::getMinimumFreeEnergy(char *sequence, char *structure) {
  int n = (int) strlen(sequence);

  /* do the linear pf fold and fill all matrices  */
  pf_linear(sequence, structure);

  /* ensemble free energy in Kcal/mol              */
  return q[iindx[1] - n];
}

/****
 * \brief TODO
 */
void MC::getProbVector(std::vector<double> &vec, size_t seq_length) {
  double tmp;

  for (size_t i = 1; i < seq_length; i++) {
    tmp = 0;
    for (size_t j = 1; j <= i - 1; j++) {
      tmp += pr[iindx[j] - i];
    }
    for (size_t j = i + 1; j <= seq_length; j++) {
      tmp += pr[iindx[i] - j];
    }
    vec.push_back(tmp);
  }
  tmp = 0;
  for (size_t j = 1; j <= seq_length - 1; j++) {
    tmp += pr[iindx[j] - seq_length];
  }
  vec.push_back(tmp);
}

/**
 * \brief TODO
 */
void MC::pf_linear(char *sequence, char *structure) {

  int n, i, j, k, l, ij, u, u1, u2, d, ii, type, type_2, tt;
  double temp, Qmax = 0;
  double qbt1, *tmp;

  double max_real;

  max_real = (sizeof(double) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  n = (int) strlen(sequence);
  //if ( n > init_length) init_pf_fold(n);  /* (re)allocate space */
  //if ((init_temp - temperature)>1e-6) update_pf_params(n);

  S = (short *) xrealloc(S, sizeof(short) * (n + 2));
  S1 = (short *) xrealloc(S1, sizeof(short) * (n + 2));
  S[0] = n;
  for (l = 1; l <= n; l++) {
    S[l] = (short) encode_char(toupper(sequence[l - 1]));
    S1[l] = alias[S[l]];
  }
  make_ptypes(S, structure);

  /* add first base at position n+1 and n'th base at position 0 */
  S[n + 1] = S[1];
  S1[n + 1] = S1[1];
  S1[0] = S1[n];

  /*array initialization ; qb,qm,q
   qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (d = 0; d <= TURN; d++)
    for (i = 1; i <= n - d; i++) {
      j = i + d;
      ij = iindx[i] - j;
      q[ij] = 1.0 * scale[d + 1];
      qb[ij] = qm[ij] = 0.0;
    }

  for (i = 1; i <= n; i++)
    qq[i] = qq1[i] = qqm[i] = qqm1[i] = 0;

  for (j = TURN + 2; j <= n; j++) {
    for (i = j - TURN - 1; i >= 1; i--) {
      /* construction of partition function of segment i,j*/
      /*firstly that given i bound to j : qb(i,j) */
      u = j - i - 1;
      ij = iindx[i] - j;
      type = ptype[ij];
      if (type != 0) {
        /*hairpin contribution*/
        qbt1 = mqtype[iindx[i + 1] - j + 1]
            * expHairpinEnergy(u, type, S1[i + 1], S1[j - 1], sequence + i - 1)
            * scale[u + 2];/* add scale[u+2] */
//				qbt1 = expHairpinEnergy(u, type, S1[i+1], S1[j-1], sequence+i-1)*scale[u+2];/* add scale[u+2] */
        /* interior loops with interior pair k,l */
        for (k = i + 1; k <= MIN(i+MAXLOOP+1,j-TURN-2); k++) {
          u1 = k - i - 1;
          for (l = MAX(k+TURN+1,j-1-MAXLOOP+u1); l < j; l++) {
            type_2 = ptype[iindx[k] - l];
            if (type_2) {
              type_2 = rtype[type_2];
              /* add *scale[u1+u2+2] */
//							qbt1 += qb[iindx[k]-l] * (scale[u1+j-l+1] *
//								expLoopEnergy(u1, j-l-1, type, type_2,
//								S1[i+1], S1[j-1], S1[k-1], S1[l+1]));
              u2 = j - l - 1;
              if (((u2 <= 1) || (mqtype[iindx[l + 1] - j + 1] != 0))
                  && ((u1 <= 1) || (mqtype[iindx[i + 1] - k + 1] != 0))
                  && ((u2 != 1) || (qtype[j - 1] != 0))
                  && ((u1 != 1) || (qtype[i + 1] != 0))) {
                //qbt1 += qb[iindx[k]-l] * (scale[u1+j-l+1] * RNAUtils::expLoopTable[u1][u2][type][type_2][S1[i+1]][S1[k-1]][S1[l+1]][S1[j-1]]);
                qbt1 += qb[iindx[k] - l]
                    * (scale[u1 + j - l + 1]
                        * expLoopEnergy(
                            u1, u2, type, type_2, S1[i + 1], S1[k - 1],
                            S1[l + 1], S1[j - 1]));
              }
            }
          }
        }
        /*multiple stem loop contribution*/
        ii = iindx[i + 1]; /* ii-k=[i+1,k-1] */
        temp = 0.0;
        for (k = i + 2; k <= j - 1; k++)
          temp += qm[ii - (k - 1)] * qqm1[k];
        tt = rtype[type];
        qbt1 += temp * expMLclosing * expMLintern[tt] * scale[2]
            * expdangle3[tt][S1[i + 1]] * expdangle5[tt][S1[j - 1]];

        qb[ij] = qbt1;
      } /* end if (type!=0) */
      else
        qb[ij] = 0.0;

      /* construction of qqm matrix containing final stem
       contributions to multiple loop partition function
       from segment i,j */
      qqm[i] = qqm1[i] * expMLbase[1];
      if (type) {
        qbt1 = qb[ij] * expMLintern[type];
        if (i > 1)
          qbt1 *= (expdangle5[type][S1[i - 1]]);
        if (j < n)
          qbt1 *= (expdangle3[type][S1[j + 1]]);
        else if (type > 2)
          qbt1 *= expTermAU;
        qqm[i] += qbt1;
      }

      /*construction of qm matrix containing multiple loop
       partition function contributions from segment i,j */
      temp = 0.0;
      ii = iindx[i]; /* ii-k=[i,k-1] */
      for (k = i + 1; k <= j; k++)
        temp += (qm[ii - (k - 1)] + expMLbase[k - i]) * qqm[k];
      qm[ij] = (temp + qqm[i]);

      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if (type) {
        if (i > 1)
          qbt1 *= (expdangle5[type][S1[i - 1]]);
        if (j < n)
          qbt1 *= (expdangle3[type][S1[j + 1]]);
        else if (type > 2)
          qbt1 *= expTermAU;
      }
      qq[i] = qq1[i] * scale[1] + qbt1;

      /*construction of partition function for segment i,j */
      temp = mqtype[iindx[i] - j] * 1.0 * scale[1 + j - i] + qq[i];
//			temp = 1.0*scale[1+j-i] + qq[i];
      for (k = i; k <= j - 1; k++)
        temp += q[ii - k] * qq[k + 1];
      q[ij] = temp;

      if (temp > Qmax) {
        Qmax = temp;
        if (Qmax > max_real / 10.)
          fprintf(stderr, "Q close to overflow: %d %d %g\n", i, j, temp);
      }
      if (temp >= max_real) {
        char msg[128];
        sprintf(msg, "overflow in pf_fold while calculating q[%d,%d]\n"
            "use larger pf_scale", i, j);
        nrerror(msg);
      }
    }
    tmp = qq1;
    qq1 = qq;
    qq = tmp;
    tmp = qqm1;
    qqm1 = qqm;
    qqm = tmp;
  }
}

/**
 * \brief Calculate base pairing probabilities
 */
void MC::pf_create_bppm(char *sequence, char *structure) {
  int n, i, j, k, l, ij, kl, ii, ll, u1, u2, type, type_2, tt, ov = 0;
  double temp, Qmax = 0, prm_MLb;
  double prmt, prmt1;
  double *tmp;

  double max_real;

  max_real = (sizeof(double) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if ((S != NULL) && (S1 != NULL)) {
    n = S[0];
    Qmax = 0;

    for (k = 1; k <= n; k++) {
      q1k[k] = q[iindx[1] - k];
      qln[k] = q[iindx[k] - n];
    }
    q1k[0] = 1.0;
    qln[n + 1] = 1.0;

    pr = q; /* recycling */

    /* 1. exterior pair i,j and initialization of pr array */
    for (i = 1; i <= n; i++) {
      for (j = i; j <= MIN(i+TURN,n); j++)
        pr[iindx[i] - j] = 0;
      for (j = i + TURN + 1; j <= n; j++) {
        ij = iindx[i] - j;
        type = ptype[ij];
        if (type && (qb[ij] > 0.)) {
          pr[ij] = q1k[i - 1] * qln[j + 1] / q1k[n];
          if (i > 1)
            pr[ij] *= expdangle5[type][S1[i - 1]];
          if (j < n)
            pr[ij] *= expdangle3[type][S1[j + 1]];
          else if (type > 2)
            pr[ij] *= expTermAU;
        } else
          pr[ij] = 0;
      }
    }

    for (l = n; l > TURN + 1; l--) {

      /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
      for (k = 1; k < l - TURN; k++) {
        kl = iindx[k] - l;
        type_2 = ptype[kl];
        type_2 = rtype[type_2];
        if (qb[kl] == 0)
          continue;

        for (i = MAX(1,k-MAXLOOP-1); i <= k - 1; i++)
          for (j = l + 1; j <= MIN(l+ MAXLOOP -k+i+2,n); j++) {
            ij = iindx[i] - j;
            type = ptype[ij];
            u1 = k - i - 1;
            u2 = j - l - 1;
            if ((pr[ij] > 0)) {
              if (((u2 <= 1) || (mqtype[iindx[l + 1] - j + 1] != 0))
                  && ((u1 <= 1) || (mqtype[iindx[i + 1] - k + 1] != 0))
                  && ((u2 != 1) || (qtype[j - 1] != 0))
                  && ((u1 != 1) || (qtype[i + 1] != 0))) {
                //pr[kl] += pr[ij] * (scale[k-i+j-l] * RNAUtils::expLoopTable[u1][u2][type][type_2][S1[i+1]][S1[k-1]][S1[l+1]][S1[j-1]]);
                pr[kl] += pr[ij]
                    * (scale[k - i + j - l]
                        * expLoopEnergy(
                            u1, u2, type, type_2, S1[i + 1], S1[k - 1],
                            S1[l + 1], S1[j - 1]));
              }
            }
          }
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
      prm_MLb = 0.;
      if (l < n)
        for (k = 2; k < l - TURN; k++) {
          i = k - 1;
          prmt = prmt1 = 0.0;

          ii = iindx[i]; /* ii-j=[i,j]     */
          ll = iindx[l + 1]; /* ll-j=[l+1,j-1] */
          tt = ptype[ii - (l + 1)];
          tt = rtype[tt];
          prmt1 = pr[ii - (l + 1)] * expMLclosing * expMLintern[tt]
              * expdangle3[tt][S1[i + 1]] * expdangle5[tt][S1[l]];
          for (j = l + 2; j <= n; j++) {
            tt = ptype[ii - j];
            tt = rtype[tt];
            prmt += pr[ii - j] * expdangle3[tt][S1[i + 1]]
                * expdangle5[tt][S1[j - 1]] * qm[ll - (j - 1)];
          }
          kl = iindx[k] - l;
          tt = ptype[kl];
          prmt *= expMLclosing * expMLintern[tt];
          prml[i] = prmt;
          prm_l[i] = prm_l1[i] * expMLbase[1] + prmt1;

          prm_MLb = prm_MLb * expMLbase[1] + prml[i];
          /* same as:    prm_MLb = 0;
           for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

          prml[i] = prml[i] + prm_l[i];

          if (qb[kl] == 0.)
            continue;

          temp = prm_MLb;

          for (i = 1; i <= k - 2; i++)
            temp += prml[i] * qm[iindx[i + 1] - (k - 1)];

          temp *= expMLintern[tt] * scale[2];
          if (k > 1)
            temp *= expdangle5[tt][S1[k - 1]];
          if (l < n)
            temp *= expdangle3[tt][S1[l + 1]];
          pr[kl] += temp;

          if (pr[kl] > Qmax) {
            Qmax = pr[kl];
            if (Qmax > max_real / 10.)
              fprintf(
                  stderr, "P close to overflow: %d %d %g %g\n", i, j, pr[kl],
                  qb[kl]);
          }
          if (pr[kl] >= max_real) {
            ov++;
            pr[kl] = FLT_MAX;
          }

        } /* end for (k=..) */
      tmp = prm_l1;
      prm_l1 = prm_l;
      prm_l = tmp;

    } /* end for (l=..)   */

    for (i = 1; i <= n; i++)
      for (j = i + TURN + 1; j <= n; j++) {
        ij = iindx[i] - j;
        pr[ij] *= qb[ij];
      }

    if (structure != NULL)
      sprintf_bppm(n, structure);
    if (ov > 0)
      fprintf(stderr, "%d overflows occurred while backtracking;\n"
          "you might try a smaller pf_scale than %g\n", ov, pf_scale);
  }/* end if((S != NULL) && (S1 != NULL))  */
  else
    nrerror(
        "bppm calculations have to be done after calling forward recursion\n");
  return;
}

/**
 * \brief TODO
 */
void MC::scale_pf_params(unsigned int length, float energy) {
  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int i, j, k, l;
  double kT, TT;
  double GT;

  kT = (temperature + K0) * GASCONST; /* kT in cal/mol  */
  TT = (temperature + K0) / (Tmeasure);

  /* scaling factors (to avoid overflows) */
  if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
    pf_scale = exp(-(-185 + (temperature - 37.) * 7.27) / kT);
    if (pf_scale < 1)
      pf_scale = 1;
  }

  double sfact = 1.07;
  double RT = (temperature + K0) * GASCONST / 1000;

  if (length > 1000)
    pf_scale = exp(-(sfact * (energy)) / RT / length);

  scale[0] = 1.;
  scale[1] = 1. / pf_scale;
  for (i = 2; i <= length; i++) {
    scale[i] = scale[i / 2] * scale[i - (i / 2)];
  }

  /* loop energies: hairpins, bulges, interior, mulit-loops */
  for (i = 0; i <= MIN(30,length); i++) {
    GT = hairpin37[i] * TT;
    exphairpin[i] = exp(-GT * 10. / kT);
  }
  for (i = 0; i <= MIN(30, MAXLOOP); i++) {
    GT = bulge37[i] * TT;
    expbulge[i] = exp(-GT * 10. / kT);
    GT = internal_loop37[i] * TT;
    expinternal[i] = exp(-GT * 10. / kT);
  }
  /* special case of size 2 interior loops (single mismatch) */
  expinternal[2] = exp(-80 * 10 / kT);

  lxc = lxc37 * TT;
  for (i = 31; i < length; i++) {
    GT = hairpin37[30] * TT + (lxc * log(i / 30.));
    exphairpin[i] = exp(-GT * 10. / kT);
  }
  for (i = 31; i <= MAXLOOP; i++) {
    GT = bulge37[30] * TT + (lxc * log(i / 30.));
    expbulge[i] = exp(-GT * 10. / kT);
    GT = internal_loop37[30] * TT + (lxc * log(i / 30.));
    expinternal[i] = exp(-GT * 10. / kT);
  }

  for (i = 0; i < 5; i++) {
    GT = F_ninio37[i] * TT;
    for (j = 0; j <= MAXLOOP; j++)
      expninio[i][j] = exp(-MIN(MAX_NINIO,j*GT) * 10 / kT);
  }
  for (i = 0; (i * 7) < strlen(Tetraloops); i++) {
    GT = TETRA_ENTH37 - (TETRA_ENTH37 - TETRA_ENERGY37[i]) * TT;
    exptetra[i] = exp(-GT * 10. / kT);
  }
  for (i = 0; (i * 5) < strlen(Triloops); i++)
    expTriloop[i] = exp(-Triloop_E37[i] * 10 / kT);

  GT = ML_closing37 * TT;
  expMLclosing = exp(-GT * 10 / kT);

  for (i = 0; i <= NBPAIRS; i++) { /* includes AU penalty */
    GT = ML_intern37 * TT;
    /* if (i>2) GT += TerminalAU; */
    expMLintern[i] = exp(-GT * 10. / kT);
  }
  expTermAU = exp(-TerminalAU * 10 / kT);

  GT = ML_BASE37 * TT;
  for (i = 0; i < length; i++) {
    expMLbase[i] = exp(-10. * i * GT / kT) * scale[i];
  }

  /* if dangles==0 just set their energy to 0,
   don't let dangle energies become > 0 (at large temps),
   but make sure go smoothly to 0                        */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= 4; j++) {
      GT = dangle5_H[i][j] - (dangle5_H[i][j] - dangle5_37[i][j]) * TT;
      expdangle5[i][j] = exp(SMOOTH(-GT) * 10. / kT);
      GT = dangle3_H[i][j] - (dangle3_H[i][j] - dangle3_37[i][j]) * TT;
      expdangle3[i][j] = exp(SMOOTH(-GT) * 10. / kT);
      if (i > 2) /* add TermAU penalty into dangle3 */
        expdangle3[i][j] *= expTermAU;
    }

  /* stacking energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++) {
      GT = enthalpies[i][j] - (enthalpies[i][j] - stack37[i][j]) * TT;
      expstack[i][j] = exp(-GT * 10 / kT);
    }

  /* mismatch energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
      for (k = 0; k < 5; k++) {
        GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchI37[i][j][k]) * TT;
        expmismatchI[i][j][k] = exp(-GT * 10.0 / kT);
        GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchH37[i][j][k]) * TT;
        expmismatchH[i][j][k] = exp(-GT * 10.0 / kT);
        GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchM37[i][j][k]) * TT;
        expmismatchM[i][j][k] = exp(-GT * 10.0 / kT);
      }

  /* interior lops of length 2 */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          GT = int11_H[i][j][k][l]
              - (int11_H[i][j][k][l] - int11_37[i][j][k][l]) * TT;
          expint11[i][j][k][l] = exp(-GT * 10. / kT);
        }
  /* interior 2x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++) {
            GT = int21_H[i][j][k][l][m]
                - (int21_H[i][j][k][l][m] - int21_37[i][j][k][l][m]) * TT;
            expint21[i][j][k][l][m] = exp(-GT * 10. / kT);
          }
        }
  /* interior 2x2 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++) {
              GT = int22_H[i][j][k][l][m][n]
                  - (int22_H[i][j][k][l][m][n] - int22_37[i][j][k][l][m][n])
                      * TT;
              expint22[i][j][k][l][m][n] = exp(-GT * 10. / kT);
            }
        }
}

/**
 * \brief TODO
 */
double MC::expHairpinEnergy(int u, int type, short si1, short sj1,
    const char *string) {
  /* compute Boltzmann weight of a hairpin loop, multiply by scale[u+2] */
  double q;
  q = exphairpin[u];
  if (u == 4) {
    char tl[7] = { 0 }, *ts;
    strncpy(tl, string, 6);
    if ((ts = strstr(Tetraloops, tl)))
      q *= exptetra[(ts - Tetraloops) / 7];
  }
  if (u == 3) {
    char tl[6] = { 0 }, *ts;
    strncpy(tl, string, 5);
    if ((ts = strstr(Triloops, tl)))
      q *= expTriloop[(ts - Triloops) / 6];
    if (type > 2)
      q *= expTermAU;
  } else
    /* no mismatches for tri-loops */
    q *= expmismatchH[type][si1][sj1];

  return q;
}

/**
 * \brief TODO
 */
double MC::expLoopEnergy(int u1, int u2, int type, int type2, short a1,
    short a2, short a3, short a4) {
  /* compute Boltzmann weight of interior loop,
   multiply by scale[u1+u2+2] for scaling */
  double z = 0;

  if ((u1 == 0) && (u2 == 0)) /* stack */
    z = expstack[type][type2];
  else {
    if ((u1 == 0) || (u2 == 0)) { /* bulge */
      int u;
      u = (u1 == 0) ? u2 : u1;
      //z = expbulge[u];
      if (u2 + u1 == 1) {
        if (u1 == 0)
          z = expbulge[u] * expstack[type][type2];
        else
          z = expbulge[u] * expstack[type][type2];
      } else {
        if (type > 2) {
          if (u1 == 0)
            z = expbulge[u] * expTermAU;
          else
            z = expbulge[u] * expTermAU;
        } else if (type2 > 2) {
          if (u1 == 0)
            z = expbulge[u] * expTermAU;
          else
            z = expbulge[u] * expTermAU;
        } else
          z = expbulge[u];
      }
    } else { /* interior loop */
      if (u1 + u2 == 2) /* size 2 is special */
        z = expint11[type][type2][a1][a4];
      else if ((u1 == 1) && (u2 == 2))
        z = expint21[type][type2][a1][a3][a4];
      else if ((u1 == 2) && (u2 == 1))
        z = expint21[type2][type][a3][a1][a2];
      else if ((u1 == 2) && (u2 == 2))
        z = expint22[type][type2][a1][a2][a3][a4];
      else {
        z = expinternal[u1 + u2] * expmismatchI[type][a1][a4]
            * expmismatchI[type2][a3][a2];
        z *= expninio[2][abs(u1 - u2)];
      }
    }
  }
  return z;
}

/**
 * \brief TODO
 */
void MC::get_arrays(unsigned int length) {
  unsigned int size, i;

  size = sizeof(double) * ((length + 1) * (length + 2) / 2);
  q = (double *) space(size);
  qb = (double *) space(size);
  qm = (double *) space(size);
  ptype = (char *) space(sizeof(char) * ((length + 1) * (length + 2) / 2));
  mqtype = (int *) space(sizeof(int) * ((length + 1) * (length + 2) / 2));

  qtype = (int *) space(sizeof(int) * (length + 1));
  q1k = (double *) space(sizeof(double) * (length + 1));
  qln = (double *) space(sizeof(double) * (length + 2));
  qq = (double *) space(sizeof(double) * (length + 2));
  qq1 = (double *) space(sizeof(double) * (length + 2));
  qqm = (double *) space(sizeof(double) * (length + 2));
  qqm1 = (double *) space(sizeof(double) * (length + 2));
  prm_l = (double *) space(sizeof(double) * (length + 2));
  prm_l1 = (double *) space(sizeof(double) * (length + 2));
  prml = (double *) space(sizeof(double) * (length + 2));
  exphairpin = (double *) space(sizeof(double) * (length + 1));
  expMLbase = (double *) space(sizeof(double) * (length + 1));
  scale = (double *) space(sizeof(double) * (length + 1));
  iindx = (int *) space(sizeof(int) * (length + 1));
  jindx = (int *) space(sizeof(int) * (length + 1));
  for (i = 1; i <= length; i++) {
    iindx[i] = ((length + 1 - i) * (length - i)) / 2 + length + 1;
    jindx[i] = (i * (i - 1)) / 2;
  }
  for (i = 0; i <= length; i++)
    qtype[i] = 1;
  for (i = 0; i < ((length + 1) * (length + 2) / 2); i++)
    mqtype[i] = 1;

}

/**
 * \brief TODO
 */
void MC::init_pf_fold(int length, float energy) {
  if (length < 1)
    nrerror("init_pf_fold: length must be greater 0");
  make_pair_matrix();
  get_arrays((unsigned) length);
  scale_pf_params((unsigned) length, energy);
}

/**
 * \brief TODO
 */
void MC::free_pf_arrays(void) {
  free(q);
  q = pr = NULL;
  free(qb);
  qb = NULL;
  free(qm);
  free(ptype);
  free(qtype);
  free(mqtype);
  free(qq);
  free(qq1);
  free(qqm);
  free(qqm1);
  free(q1k);
  free(qln);
  free(prm_l);
  free(prm_l1);
  free(prml);
  free(exphairpin);
  free(expMLbase);
  free(scale);
  free(iindx);
  free(jindx);
  free(S);
  S = NULL;
  free(S1);
  S1 = NULL;
}

/**
 * \brief TODO
 */
void MC::update_pf_params(int length, float energy) {
  make_pair_matrix();
  scale_pf_params((unsigned) length, energy);
}

/**
 * \brief TODO
 */
char MC::bppm_symbol(float *x) {

  if (x[0] > 0.667)
    return '.';
  if (x[1] > 0.667)
    return '(';
  if (x[2] > 0.667)
    return ')';
  if ((x[1] + x[2]) > x[0]) {
    if ((x[1] / (x[1] + x[2])) > 0.667)
      return '{';
    if ((x[2] / (x[1] + x[2])) > 0.667)
      return '}';
    else
      return '|';
  }
  if (x[0] > (x[1] + x[2]))
    return ',';
  return ':';
}

/**
 * \brief TODO
 */
void MC::sprintf_bppm(int length, char *structure) {
  int i, j;
  float P[3]; /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */

  for (j = 1; j <= length; j++) {
    P[0] = 1.0;
    P[1] = P[2] = 0.0;
    for (i = 1; i < j; i++) {
      P[2] += pr[iindx[i] - j]; /* j is paired downstream */
      P[0] -= pr[iindx[i] - j]; /* j is unpaired */
    }
    for (i = j + 1; i <= length; i++) {
      P[1] += pr[iindx[j] - i]; /* j is paired upstream */
      P[0] -= pr[iindx[j] - i]; /* j is unpaired */
    }
    structure[j - 1] = bppm_symbol(P);
  }
  structure[length] = '\0';
}

/**
 * \brief TODO
 */
void MC::make_ptypes(const short *S, const char *structure) {
  int n, i, j, k, l;

  n = S[0];
  for (k = 1; k < n - TURN; k++)
    for (l = 1; l <= 2; l++) {
      int type, ntype = 0;
      i = k;
      j = i + TURN + l;
      if (j > n)
        continue;
      type = pair[S[i]][S[j]];
      while ((i >= 1) && (j <= n)) {
        if ((i > 1) && (j < n))
          ntype = pair[S[i - 1]][S[j + 1]];
        qb[iindx[i] - j] = 0.;
        ptype[iindx[i] - j] = (char) type;
        type = ntype;
        i--;
        j++;
      }
    }

  if (structure != NULL) {
    int hx, *stack;
    char type;
    stack = (int *) space(sizeof(int) * (n + 1));

    for (hx = 0, j = 1; j <= n; j++) {
      switch (structure[j - 1]) {
      case 'x': /* can't pair */
        for (l = 1; l < j - TURN; l++)
          ptype[iindx[l] - j] = 0;
        for (l = j + TURN + 1; l <= n; l++)
          ptype[iindx[j] - l] = 0;
        break;
      case '(':
        stack[hx++] = j;
        /* fallthrough */
      case '<': /* pairs upstream */
        for (l = 1; l < j - TURN; l++)
          ptype[iindx[l] - j] = 0;
        break;
      case ')':
        if (hx <= 0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack[--hx];
        type = ptype[iindx[i] - j];
        /* don't allow pairs i<k<j<l */
        for (k = i; k <= j; k++)
          for (l = j; l <= n; l++)
            ptype[iindx[k] - l] = 0;
        /* don't allow pairs k<i<l<j */
        for (k = 1; k <= i; k++)
          for (l = i; l <= j; l++)
            ptype[iindx[k] - l] = 0;
        ptype[iindx[i] - j] = (type == 0) ? 7 : type;
        /* fallthrough */
      case '>': /* pairs downstream */
        for (l = j + TURN + 1; l <= n; l++)
          ptype[iindx[j] - l] = 0;
        break;
      case '^': /* must pair with something */
        qtype[j] = 0;
        break;
      }
    }
    if (hx != 0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in constraint string");
    }
    free(stack);
  }
  for (i = 1; i <= n; i++)
    mqtype[iindx[i] - i] = qtype[i];
  for (i = 1; i < n; i++) {
    for (j = i + 1; j <= n; j++)
      mqtype[iindx[i] - j] = mqtype[iindx[i] - j + 1] * qtype[j];
  }
}

/**
 * \brief TODO
 */
double MC::getTemperature() {
  return temperature;
}

/****
 * \brief write the base-pair probability matrix(?) from this MC object to
 *        the given file.
 * \todo replace with iostreams
 */
void MC::printMatrix(string fname, int n, int seql, int pos) {
  FILE *fout;
  fout = fopen(fname.c_str(), "a+");
  double tmp;
  int i, j, ij;

  for (i = 1; i < seql; i++) {
    tmp = 0;
    for (j = 1; j <= i - 1; j++) {
      ij = iindx[j] - i;
      tmp += pr[ij];
    }
    for (j = i + 1; j <= seql; j++) {
      ij = iindx[i] - j;
      tmp += pr[ij];
    }
    fprintf(fout, "%g,", tmp);
  }
  tmp = 0;
  for (j = 1; j <= seql - 1; j++) {
    ij = iindx[j] - seql;
    tmp += pr[ij];
  }
  fprintf(fout, "%g\n", tmp);

  fclose(fout);
}

