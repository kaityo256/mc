//----------------------------------------------------------------------
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
//----------------------------------------------------------------------
const int L = 64;       // System Size
const int N = L*L;      // Number of spins
const int T_LOOP = 1000; // Thermalization loop
const int S_LOOP = 1000; // Number of average
const int O_LOOP = 10;  // Number of observation at a temperature
const int ND = 20;      // Number of temperatures to be observed
std::vector<int> spin(N), cluster(N), flip(N);
std::mt19937 mt(1);
std::uniform_real_distribution<double> ud(0.0, 1.0);
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
connect_bond(int ix1, int iy1, int ix2, int iy2, const double p) {
  if (ix2 == L)ix2 = 0;
  if (iy2 == L)iy2 = 0;
  if (spin[ix1 + iy1 * L] * spin[ix2 + iy2 * L] < 0)return;
  if (ud(mt) >= p)return;
  const int i1 = ix1 + iy1 * L;
  const int i2 = ix2 + iy2 * L;
  const int c1 = get_cluster_number(i1);
  const int c2 = get_cluster_number(i2);
  if (c1 < c2) {
    cluster[c2] = c1;
  } else {
    cluster[c1] = c2;
  }
}
//----------------------------------------------------------------------
void
cluster_flip(double beta) {
  for (int i = 0; i < N; i++) {
    cluster[i] = i;
    flip[i] = static_cast<int>(ud(mt) * 2.0) * 2 - 1;
  }
  const double p = 1.0 - exp(-2.0 * beta);
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      connect_bond(ix, iy, ix + 1, iy, p );
      connect_bond(ix, iy, ix, iy + 1, p );
    }
  }
  for (int i=0;i < N;i++){
    const int c = get_cluster_number(i);
    spin[i] = flip[c];
  }
}
//----------------------------------------------------------------------
void
single_flip(double beta) {
  for(int i=0;i<N;i++){
    int ix = i % L;
    int iy = i /L;
    double de = 0.0;
    if (ix != 0)de += spin[(ix - 1) + iy * L];
    else de += spin[(L - 1) + iy * L];
    if (iy != 0)de += spin[ix + (iy - 1)*L];
    else de += spin[ix + (L - 1)*L];
    if (ix != L - 1)de += spin[(ix + 1) + iy*L];
    else de += spin[0 + iy * L];
    if (iy != L - 1)de += spin[ix + (iy + 1)*L];
    else de += spin[ix + 0 * L];
    de = de * 2.0 * spin[ix + iy * L];
    if (de < 0 || exp(-beta * de) > ud(mt)) {
      spin[i] = - spin[i];
    }
  }
}
//----------------------------------------------------------------------
void
mc_onestep(double beta) {
#ifdef CLUSTER
  cluster_flip(beta);
#else
  single_flip(beta);
#endif
}
//----------------------------------------------------------------------
// Observation
//----------------------------------------------------------------------
double
magnetization_square(void) {
  double m = std::accumulate(spin.begin(),spin.end(),0.0);
  m /= static_cast<double>(N);
  return m * m;
}
//----------------------------------------------------------------------
double
energy(void) {
  double e = 0.0;
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      const int s = spin[ix + iy* L];
      if (ix != L - 1) e -= s * spin[ix + 1 + iy*L];
      else        e -= s * spin[0 + iy* L];
      if (iy != L - 1) e -= s * spin[ix +(iy + 1)*L];
      else        e -= s * spin[ix + 0 * L];
    }
  }
  e /= static_cast<double>(N);
  return e;
}
//----------------------------------------------------------------------
void
domc(void) {
  double bs = 0.3; //Start inverse temperature  (The highest temperature)
  double be = 0.6; //End inverse temperature (The lowest temperature)
  for (int i = 0; i < ND; i++) {
    double beta = bs + (be - bs) * (double)i / (double)ND;
    std::cerr << "beta = " << beta << std::endl;
    std::fill(spin.begin(), spin.end(), 1);
    for (int j = 0; j < T_LOOP; j++) {
      mc_onestep(beta);
    }
    for (int j = 0; j < O_LOOP; j++) {
      double m2 = 0.0;
      double e_sum = 0.0;
      double e2_sum = 0.0;
      for (int k = 0; k < S_LOOP; k++) {
        mc_onestep(beta);
        m2 += magnetization_square();
        const double e = energy();
        e_sum += e;
        e2_sum += (e * e);
      }
      m2 /= static_cast<double>(S_LOOP);
      e_sum /= static_cast<double>(S_LOOP);
      e2_sum /= static_cast<double>(S_LOOP);
      const double c = (e2_sum - e_sum * e_sum) * beta * beta * N;
      std::cout << beta << " " << m2 << " " << c << std::endl;
    }
  }
}
//----------------------------------------------------------------------
int
main(void) {
#ifdef CLUSTER
  std::cerr << "Swendsen-Wang" << std::endl;
#else
  std::cerr << "Single Flip" << std::endl;
#endif
  std::cout << "# beta, Magnetization^2, Specific Heat" << std::endl;
  domc();
}
//----------------------------------------------------------------------
