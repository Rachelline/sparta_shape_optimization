/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "spatype.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math_extra.h"
#include <type_traits>

using namespace SPARTA_NS;

namespace MathExtra {

namespace detail {

  // trivially-copyable types (e.g. double): memcpy/memmove is safe and fast
  template <unsigned N, typename T>
  inline typename std::enable_if<std::is_trivially_copyable<T>::value>::type
  swap_row(T *a, T *b)
  {
    T tempv[N];
    memcpy(tempv,a,N*sizeof(T));
    memmove(a,b,N*sizeof(T));
    memcpy(b,tempv,N*sizeof(T));
  }

  // non-trivially-copyable types (e.g. AD active types): must go through
  // assignment so constructors/tape bookkeeping run correctly
  template <unsigned N, typename T>
  inline typename std::enable_if<!std::is_trivially_copyable<T>::value>::type
  swap_row(T *a, T *b)
  {
    for (unsigned k = 0; k < N; k++) {
      T tempv = a[k];
      a[k] = b[k];
      b[k] = tempv;
    }
  }

}

/* ----------------------------------------------------------------------
   solve Ax = b or M ans = v
   use gaussian elimination & partial pivoting on matrix
   copied from lammps
------------------------------------------------------------------------- */

template <typename T>
int mldivide3(const T m[3][3], const T *v, T *ans)
{
  // create augmented matrix for pivoting

  T aug[3][4];
  for (unsigned i = 0; i < 3; i++) {
    aug[i][3] = v[i];
    for (unsigned j = 0; j < 3; j++) aug[i][j] = m[i][j];
  }

  for (unsigned i = 0; i < 2; i++) {
    unsigned p = i;
    for (unsigned j = i+1; j < 3; j++) {
      if (fabs(aug[j][i]) > fabs(aug[i][i])) {
        detail::swap_row<4>(aug[i],aug[j]);
      }
    }

    while (p < 3 && aug[p][i] == 0.0) p++;

    if (p == 3) return 1;
    else
      if (p != i) {
        detail::swap_row<4>(aug[i],aug[p]);
      }

    for (unsigned j = i+1; j < 3; j++) {
      T n = aug[j][i]/aug[i][i];
      for (unsigned k=i+1; k<4; k++) aug[j][k]-=n*aug[i][k];
    }
  }

  if (aug[2][2] == 0.0) return 1;

  // back substitution

  ans[2] = aug[2][3]/aug[2][2];
  for (int i = 1; i >= 0; i--) {
    T sumax = 0.0;
    for (unsigned j = i+1; j < 3; j++) sumax += aug[i][j]*ans[j];
    ans[i] = (aug[i][3]-sumax) / aug[i][i];
  }

  return 0;
}


/* ----------------------------------------------------------------------
   solve Ax = b or M ans = v
   use gaussian elimination & partial pivoting on matrix
   returns 1 if fails? and 0 if succeeds?
   copied from LAMMPS
------------------------------------------------------------------------- */

template <typename T>
int mldivide4(const T m[4][4], const T *v, T *ans)
{
  // create augmented matrix for pivoting

  T aug[4][5];
  for (unsigned i = 0; i < 4; i++) {
    aug[i][4] = v[i];
    for (unsigned j = 0; j < 4; j++) aug[i][j] = m[i][j];
  }

  for (unsigned i = 0; i < 3; i++) {
    unsigned p = i;

    // swaps rows so largest ones are at top

    for (unsigned j = i+1; j < 4; j++) {
      if (fabs(aug[j][i]) > fabs(aug[i][i])) {
        detail::swap_row<5>(aug[i],aug[j]);
      }
    }

    // checks how many zeros in column i

    while (p < 4 && aug[p][i] == 0.0) p++;

    // if all zero no solution

    if (p == 4) return 1;
    else {
      if (p != i) {
        detail::swap_row<5>(aug[i],aug[p]);
      }
    }

    // Gaussian elimination

    for (unsigned j = i+1; j < 4; j++) {
      T n = aug[j][i]/aug[i][i];
      for (unsigned k = i+1; k < 5; k++) aug[j][k]-=n*aug[i][k];
    }
  }

  // if last diagonal is zero, no solution

  if (aug[3][3] == 0.0) return 1;

  // back substitution

  ans[3] = aug[3][4]/aug[3][3];
  for (int i = 2; i >= 0; i--) {
    T sumax = 0.0;
    for (unsigned j = i+1; j < 4; j++) sumax += aug[i][j]*ans[j];
    ans[i] = (aug[i][4]-sumax) / aug[i][i];
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, Nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 1 to Nmax,
     (3) i* = i to Nmax, (4) *j = 1 to j, (5) i*j = i to j
   return nlo,nhi
   return 0 if successful
   return 1 if numeric values are out of lower/upper bounds
------------------------------------------------------------------------- */

int bounds(char *str, int nmax, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = nhi = atoi(str);
  } else if (strlen(str) == 1) {
    nlo = 1;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = 1;
    nhi = atoi(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = atoi(str);
    nhi = nmax;
  } else {
    nlo = atoi(str);
    nhi = atoi(ptr+1);
  }

  if (nlo < 1 || nhi > nmax) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   convert a 64-bit integer into a string like "9.36B"
   K = thousand, M = million, B = billion, T = trillion, P = peta, E = exa
   for easier-to-understand output
------------------------------------------------------------------------- */

char *num2str(bigint n, char *outstr)
{
  if (n < 100000) sprintf(outstr,"(%1.3gK)",1.0e-3*n);
  else if (n < 1000000000) sprintf(outstr,"(%1.3gM)",1.0e-6*n);
  else if (n < 1000000000000) sprintf(outstr,"(%1.3gB)",1.0e-9*n);
  else if (n < 1000000000000000) sprintf(outstr,"(%1.3gT)",1.0e-12*n);
  else if (n < 1000000000000000000) sprintf(outstr,"(%1.3gP)",1.0e-15*n);
  else sprintf(outstr,"(%1.3gE)",1.0e-18*n);
  return outstr;
}

// explicit instantiations

template int mldivide3<double>(const double m[3][3], const double *v, double *ans);
template int mldivide4<double>(const double m[4][4], const double *v, double *ans);

}
