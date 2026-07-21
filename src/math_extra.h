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

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#ifndef SPARTA_MATH_EXTRA_H
#define SPARTA_MATH_EXTRA_H

#include "spatype.h"
#include "math.h"
#include "stdio.h"
#include "string.h"

namespace MathExtra {

  // 3 vector operations

  template <typename T>
  void norm3(T *v);
  template <typename T>
  void normalize3(const T *v, T *ans);
  template <typename T>
  void snorm3(const T, T *v);
  template <typename T>
  void snormalize3(const T, const T *v, T *ans);
  template <typename T>
  void negate3(T *v);
  template <typename T>
  void scale3(T s, T *v);
  template <typename T>
  void scale3(T s, const T *v, T *ans);
  template <typename T>
  void axpy3(T alpha, const T *x, T *y);
  template <typename T>
  void axpy3(T alpha, const T *x, const T *y,
             T *ynew);
  template <typename T>
  void add3(const T *v1, const T *v2, T *ans);
  template <typename T>
  void sub3(const T *v1, const T *v2, T *ans);
  template <typename T>
  T len3(const T *v);
  template <typename T>
  T lensq3(const T *v);
  template <typename T>
  T dot3(const T *v1, const T *v2);
  template <typename T>
  void cross3(const T *v1, const T *v2, T *ans);
  template <typename T>
  void reflect3(T *v, const T *norm);

  // 3x3 matrix operations

  template <typename T>
  T det3(const T mat[3][3]);
  template <typename T>
  void diag_times3(const T *diagonal, const T mat[3][3],
                    T ans[3][3]);
  template <typename T>
  void plus3(const T m[3][3], const T m2[3][3],
             T ans[3][3]);
  template <typename T>
  void times3(const T m[3][3], const T m2[3][3],
              T ans[3][3]);
  template <typename T>
  void transpose_times3(const T mat1[3][3],
                         const T mat2[3][3],
                         T ans[3][3]);
  template <typename T>
  void times3_transpose(const T mat1[3][3],
                         const T mat2[3][3],
                         T ans[3][3]);
  template <typename T>
  void invert3(const T mat[3][3], T ans[3][3]);
  template <typename T>
  void matvec(const T mat[3][3], const T *vec, T *ans);
  template <typename T>
  void matvec(const T *ex, const T *ey, const T *ez,
              const T *vec, T *ans);
  template <typename T>
  void transpose_matvec(const T mat[3][3], const T *vec,
                         T *ans);
  template <typename T>
  void transpose_matvec(const T *ex, const T *ey,
                         const T *ez, const T *v,
                         T *ans);
  template <typename T>
  void transpose_diag3(const T mat[3][3], const T *vec,
                        T ans[3][3]);
  template <typename T>
  void vecmat(const T *v, const T m[3][3], T *ans);
  template <typename T>
  void scalar_times3(const T f, T m[3][3]);

  template <typename T>
  int mldivide3(const T mat[3][3], const T *vec, T *ans);
  template <typename T>
  int mldivide4(const T mat[4][4], const T *vec, T *ans);

  // quaternion operations

  template <typename T>
  void axisangle_to_quat(const T *v, const T angle,
                          T *quat);
  template <typename T>
  void quat_to_mat(const T *quat, T mat[3][3]);

  // misc methods

  int bounds(char *, int, int &, int &);
  char *num2str(SPARTA_NS::bigint, char *);
}

/* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::norm3(T *v)
{
  T scale = T(1.0)/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

/* ----------------------------------------------------------------------
   normalize a vector, return in ans
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::normalize3(const T *v, T *ans)
{
  T scale = T(1.0)/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length in place
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::snorm3(const T length, T *v)
{
  T scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::snormalize3(const T length, const T *v, T *ans)
{
  T scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   negate vector v in place
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::negate3(T *v)
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

/* ----------------------------------------------------------------------
   scale vector v by s in place
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::scale3(T s, T *v)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

/* ----------------------------------------------------------------------
   scale vector v by s, return in ans
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::scale3(T s, const T *v, T *ans)
{
  ans[0] = s*v[0];
  ans[1] = s*v[1];
  ans[2] = s*v[2];
}

/* ----------------------------------------------------------------------
   axpy: y = alpha*x + y
   y is replaced by result
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::axpy3(T alpha, const T *x, T *y)
{
  y[0] += alpha*x[0];
  y[1] += alpha*x[1];
  y[2] += alpha*x[2];
}

/* ----------------------------------------------------------------------
   axpy: ynew = alpha*x + y
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::axpy3(T alpha, const T *x, const T *y,
                      T *ynew)
{
  ynew[0] += alpha*x[0] + y[0];
  ynew[1] += alpha*x[1] + y[1];
  ynew[2] += alpha*x[2] + y[2];
}

/* ----------------------------------------------------------------------
   ans = v1 + v2
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::add3(const T *v1, const T *v2, T *ans)
{
  ans[0] = v1[0] + v2[0];
  ans[1] = v1[1] + v2[1];
  ans[2] = v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = v1 - v2
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::sub3(const T *v1, const T *v2, T *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

template <typename T>
T MathExtra::len3(const T *v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

template <typename T>
T MathExtra::lensq3(const T *v)
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

template <typename T>
T MathExtra::dot3(const T *v1, const T *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::cross3(const T *v1, const T *v2, T *ans)
{
  ans[0] = v1[1]*v2[2] - v1[2]*v2[1];
  ans[1] = v1[2]*v2[0] - v1[0]*v2[2];
  ans[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   reflect vector v around unit normal n
   return updated v of same length = v - 2(v dot n)n
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::reflect3(T *v, const T *n)
{
  T dot = dot3(v,n);
  v[0] -= 2.0*dot*n[0];
  v[1] -= 2.0*dot*n[1];
  v[2] -= 2.0*dot*n[2];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

template <typename T>
T MathExtra::det3(const T m[3][3])
{
  T ans = m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] -
    m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1] +
    m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::diag_times3(const T *d, const T m[3][3],
                            T ans[3][3])
{
  ans[0][0] = d[0]*m[0][0];
  ans[0][1] = d[0]*m[0][1];
  ans[0][2] = d[0]*m[0][2];
  ans[1][0] = d[1]*m[1][0];
  ans[1][1] = d[1]*m[1][1];
  ans[1][2] = d[1]*m[1][2];
  ans[2][0] = d[2]*m[2][0];
  ans[2][1] = d[2]*m[2][1];
  ans[2][2] = d[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::plus3(const T m[3][3], const T m2[3][3],
                      T ans[3][3])
{
  ans[0][0] = m[0][0]+m2[0][0];
  ans[0][1] = m[0][1]+m2[0][1];
  ans[0][2] = m[0][2]+m2[0][2];
  ans[1][0] = m[1][0]+m2[1][0];
  ans[1][1] = m[1][1]+m2[1][1];
  ans[1][2] = m[1][2]+m2[1][2];
  ans[2][0] = m[2][0]+m2[2][0];
  ans[2][1] = m[2][1]+m2[2][1];
  ans[2][2] = m[2][2]+m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times mat2
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::times3(const T m[3][3], const T m2[3][3],
                       T ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[1][0] + m[0][2]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[0][1]*m2[1][1] + m[0][2]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[0][1]*m2[1][2] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[1][0] + m[1][2]*m2[2][0];
  ans[1][1] = m[1][0]*m2[0][1] + m[1][1]*m2[1][1] + m[1][2]*m2[2][1];
  ans[1][2] = m[1][0]*m2[0][2] + m[1][1]*m2[1][2] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[2][0]*m2[0][1] + m[2][1]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[2][0]*m2[0][2] + m[2][1]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::transpose_times3(const T m[3][3], const T m2[3][3],
                                 T ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[1][0]*m2[1][0] + m[2][0]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[1][0]*m2[1][1] + m[2][0]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[1][0]*m2[1][2] + m[2][0]*m2[2][2];
  ans[1][0] = m[0][1]*m2[0][0] + m[1][1]*m2[1][0] + m[2][1]*m2[2][0];
  ans[1][1] = m[0][1]*m2[0][1] + m[1][1]*m2[1][1] + m[2][1]*m2[2][1];
  ans[1][2] = m[0][1]*m2[0][2] + m[1][1]*m2[1][2] + m[2][1]*m2[2][2];
  ans[2][0] = m[0][2]*m2[0][0] + m[1][2]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[0][2]*m2[0][1] + m[1][2]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[0][2]*m2[0][2] + m[1][2]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times transpose of mat2
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::times3_transpose(const T m[3][3], const T m2[3][3],
                                 T ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[0][1] + m[0][2]*m2[0][2];
  ans[0][1] = m[0][0]*m2[1][0] + m[0][1]*m2[1][1] + m[0][2]*m2[1][2];
  ans[0][2] = m[0][0]*m2[2][0] + m[0][1]*m2[2][1] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[0][1] + m[1][2]*m2[0][2];
  ans[1][1] = m[1][0]*m2[1][0] + m[1][1]*m2[1][1] + m[1][2]*m2[1][2];
  ans[1][2] = m[1][0]*m2[2][0] + m[1][1]*m2[2][1] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[0][1] + m[2][2]*m2[0][2];
  ans[2][1] = m[2][0]*m2[1][0] + m[2][1]*m2[1][1] + m[2][2]*m2[1][2];
  ans[2][2] = m[2][0]*m2[2][0] + m[2][1]*m2[2][1] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   invert a matrix
   does NOT checks for singular or badly scaled matrix
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::invert3(const T m[3][3], T ans[3][3])
{
  T den = m[0][0]*m[1][1]*m[2][2]-m[0][0]*m[1][2]*m[2][1];
  den += -m[1][0]*m[0][1]*m[2][2]+m[1][0]*m[0][2]*m[2][1];
  den += m[2][0]*m[0][1]*m[1][2]-m[2][0]*m[0][2]*m[1][1];

  ans[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1]) / den;
  ans[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1]) / den;
  ans[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1]) / den;
  ans[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0]) / den;
  ans[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0]) / den;
  ans[1][2] = -(m[0][0]*m[1][2]-m[0][2]*m[1][0]) / den;
  ans[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0]) / den;
  ans[2][1] = -(m[0][0]*m[2][1]-m[0][1]*m[2][0]) / den;
  ans[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) / den;
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::matvec(const T m[3][3], const T *v, T *ans)
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::matvec(const T *ex, const T *ey, const T *ez,
                       const T *v, T *ans)
{
  ans[0] = ex[0]*v[0] + ey[0]*v[1] + ez[0]*v[2];
  ans[1] = ex[1]*v[0] + ey[1]*v[1] + ez[1]*v[2];
  ans[2] = ex[2]*v[0] + ey[2]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::transpose_matvec(const T m[3][3], const T *v,
                                 T *ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::transpose_matvec(const T *ex, const T *ey,
                                 const T *ez, const T *v,
                                 T *ans)
{
  ans[0] = ex[0]*v[0] + ex[1]*v[1] + ex[2]*v[2];
  ans[1] = ey[0]*v[0] + ey[1]*v[1] + ey[2]*v[2];
  ans[2] = ez[0]*v[0] + ez[1]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::transpose_diag3(const T m[3][3], const T *d,
                                T ans[3][3])
{
  ans[0][0] = m[0][0]*d[0];
  ans[0][1] = m[1][0]*d[1];
  ans[0][2] = m[2][0]*d[2];
  ans[1][0] = m[0][1]*d[0];
  ans[1][1] = m[1][1]*d[1];
  ans[1][2] = m[2][1]*d[2];
  ans[2][0] = m[0][2]*d[0];
  ans[2][1] = m[1][2]*d[1];
  ans[2][2] = m[2][2]*d[2];
}

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::vecmat(const T *v, const T m[3][3], T *ans)
{
  ans[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0];
  ans[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1];
  ans[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   matrix times scalar, in place
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::scalar_times3(const T f, T m[3][3])
{
  m[0][0] *= f; m[0][1] *= f; m[0][2] *= f;
  m[1][0] *= f; m[1][1] *= f; m[1][2] *= f;
  m[2][0] *= f; m[2][1] *= f; m[2][2] *= f;
}

/* ----------------------------------------------------------------------
   compute quaternion from axis-angle rotation
   v MUST be a unit vector
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::axisangle_to_quat(const T *v, const T angle,
                                  T *quat)
{
  T halfa = 0.5*angle;
  T sina = sin(halfa);
  quat[0] = cos(halfa);
  quat[1] = v[0]*sina;
  quat[2] = v[1]*sina;
  quat[3] = v[2]*sina;
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion
   quat = [w i j k]
------------------------------------------------------------------------- */

template <typename T>
void MathExtra::quat_to_mat(const T *quat, T mat[3][3])
{
  T w2 = quat[0]*quat[0];
  T i2 = quat[1]*quat[1];
  T j2 = quat[2]*quat[2];
  T k2 = quat[3]*quat[3];
  T twoij = 2.0*quat[1]*quat[2];
  T twoik = 2.0*quat[1]*quat[3];
  T twojk = 2.0*quat[2]*quat[3];
  T twoiw = 2.0*quat[1]*quat[0];
  T twojw = 2.0*quat[2]*quat[0];
  T twokw = 2.0*quat[3]*quat[0];

  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;

  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;

  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}

#endif
