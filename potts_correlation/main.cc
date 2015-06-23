//----------------------------------------------------------------------
//          Copyright H. Watanabe 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//     http://www.boost.org/LICENSE_1_0.txt)
//----------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fftw3.h>
#include "mt.h"
//----------------------------------------------------------------------
const int LX = 128;
const int LY = 128;
const int N = LX * LY;
const int Q = 3;
int T_LOOP = 1000;
int O_LOOP = 10000;
int spin[N], cluster[N], newspin[N];
double sofk[N];
//----------------------------------------------------------------------
int
get_cluster_number(int index) {
  int i = index;
  while (i != cluster[i]) {
    i = cluster[i];
  }
  return i;
}
//----------------------------------------------------------------------
void
connect(const int i1, const int i2, double p) {
  if (spin[i1] != spin[i2]) return;
  if (MT::GetDouble() > p) return;
  const int c1 = get_cluster_number(i1);
  const int c2 = get_cluster_number(i2);
  if (c1 < c2) {
    cluster[c2] = c1;
  } else {
    cluster[c1] = c2;
  }
}
//----------------------------------------------------------------------
inline
int
pos2index(int ix, int iy) {
  ix = ix % LX;
  iy = iy % LY;
  return ix + iy * LX;
}
//----------------------------------------------------------------------
void
sw_flip(const double beta) {
  const double p = 1.0 - exp(-beta);
  for (int i = 0; i < N; i++) {
    cluster[i] = i;
    newspin[i] = static_cast<int>(MT::GetDouble() * Q);
  }
  for (int i = 0; i < N; i++) {
    const int ix = i % LX;
    const int iy = i / LX;
    connect(i, pos2index(ix + 1, iy), p);
    connect(i, pos2index(ix, iy + 1), p);
  }
  for (int i = 0; i < N; i++) {
    const int c = get_cluster_number(i);
    spin[i] = newspin[c];
  }
}
//----------------------------------------------------------------------
void
calc_correlation(fftw_complex *in, fftw_complex *out, fftw_plan &p) {
  const double q = static_cast<double>(Q);
  for (int i = 0; i < N; i++) {
    if(spin[i] == 0){
      in[i][0] = 1.0;
    }else{
      in[i][0] = -1.0/(q-1);
    }
    in[i][1] = 0.0;
  }
  fftw_execute(p);
  const double scale = 1.0 / static_cast<double>(N);
  for (int i = 0; i < N; i++) {
    const double re = out[i][0];
    const double im = out[i][1];
    const double v = (re * re + im * im) * scale;
    sofk[i] += v;
  }
}
//----------------------------------------------------------------------
int
main(void) {
  const double bc = log(1.0 + sqrt(Q));
  fftw_complex *in = NULL;
  fftw_complex *out = NULL;
  fftw_plan p = NULL;
  const size_t s = sizeof(fftw_complex) * N;
  in = (fftw_complex*) fftw_malloc(s);
  out = (fftw_complex*) fftw_malloc(s);
  p = fftw_plan_dft_2d(LX, LY, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  for (int i = 0; i < N; i++) {
    spin[i] = 0;
    sofk[i] = 0.0;
  }
  //Thermalization
  for (int i = 0; i < T_LOOP; i++) {
    sw_flip(bc);
  }
  //Observation
  for (int i = 0; i < O_LOOP; i++) {
    sw_flip(bc);
    calc_correlation(in, out, p);
  }

  for (int i = 0; i < N; i++) {
    in[i][0] = sofk[i] / static_cast<double>(O_LOOP);
    in[i][1] = 0.0;
  }
  p = fftw_plan_dft_2d(LX, LY, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  printf("# Q=%d\n", Q);
  printf("# T_LOOP=%d\n", T_LOOP);
  printf("# O_LOOP=%d\n", O_LOOP);
  const int l = std::min(LX,LY)/2;
  for (int i = 0; i < l; i++) {
    for (int j = 0; j < l; j++) {
      const int index = j * LX + i;
      const double re = out[index][0];
      const double dx = static_cast<double>(i);
      const double dy = static_cast<double>(j);
      printf("%f %f\n", sqrt(dx*dx+dy*dy), re/static_cast<double>(N)*(Q-1));
    }
  }
  fftw_free(in);
  fftw_free(out);
}
//----------------------------------------------------------------------
