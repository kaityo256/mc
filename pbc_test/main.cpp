#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

const int L = 64;    // System Size
const int N = L * L; // Number of spins
double e_table[5];   // Energy Table

std::vector<int> spin(N); //Spins

std::mt19937 mt(1);
std::uniform_real_distribution<double> ud(0.0, 1.0);

int pos2index_mod(int ix, int iy) {
  ix = (ix + L) % L;
  iy = (iy + L) % L;
  return ix + iy * L;
}

int pos2index_if(int ix, int iy) {
  if (ix < 0) ix += L;
  if (ix >= L) ix -= L;
  if (iy < 0) iy += L;
  if (iy >= L) iy -= L;
  return ix + iy * L;
}

// For Neighbor list
std::vector<std::array<int, 4>> neighbor(N);

void init_neighbors() {
  for (int i = 0; i < N; i++) {
    int ix = i % L;
    int iy = i / L;
    neighbor[i][0] = pos2index_mod(ix + 1, iy);
    neighbor[i][1] = pos2index_mod(ix - 1, iy);
    neighbor[i][2] = pos2index_mod(ix, iy + 1);
    neighbor[i][3] = pos2index_mod(ix, iy - 1);
  }
}

// Single Flip Implementation (Metropolis Algorithm)

/*
  Adjust PBC by if
  Calculate Energy Explicitly
*/
void single_flip_calc_if(double beta) {
  for (int i = 0; i < N; i++) {
    int ix = i % L;
    int iy = i / L;
    int ns = 0;
    ns += spin[pos2index_if(ix + 1, iy)];
    ns += spin[pos2index_if(ix - 1, iy)];
    ns += spin[pos2index_if(ix, iy + 1)];
    ns += spin[pos2index_if(ix, iy - 1)];
    ns *= spin[i];
    if (ns < 0 || exp(-2.0 * ns * beta) > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
}

/*
  Adjust PBC by mod
  Calculate Energy Explicitly
*/
void single_flip_calc_mod(double beta) {
  for (int i = 0; i < N; i++) {
    int ix = i % L;
    int iy = i / L;
    int ns = 0;
    ns += spin[pos2index_mod(ix + 1, iy)];
    ns += spin[pos2index_mod(ix - 1, iy)];
    ns += spin[pos2index_mod(ix, iy + 1)];
    ns += spin[pos2index_mod(ix, iy - 1)];
    ns *= spin[i];
    if (ns < 0 || exp(-2.0 * ns * beta) > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
}

/*
  Adjust PBC by mod
  Use a table for Energy
*/
void single_flip_table_mod(double beta) {
  for (int i = 0; i < N; i++) {
    int ix = i % L;
    int iy = i / L;
    int ns = 0;
    ns += spin[pos2index_mod(ix + 1, iy)];
    ns += spin[pos2index_mod(ix - 1, iy)];
    ns += spin[pos2index_mod(ix, iy + 1)];
    ns += spin[pos2index_mod(ix, iy - 1)];
    ns *= spin[i];
    if (ns < 0 || e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
}

/*
  Adjust PBC by mod
  Use a table for Energy
  omit (if ns < 0)
*/
void single_flip_table_mod_tableonly(double beta) {
  for (int i = 0; i < N; i++) {
    int ix = i % L;
    int iy = i / L;
    int ns = 0;
    ns += spin[pos2index_mod(ix + 1, iy)];
    ns += spin[pos2index_mod(ix - 1, iy)];
    ns += spin[pos2index_mod(ix, iy + 1)];
    ns += spin[pos2index_mod(ix, iy - 1)];
    ns *= spin[i];
    if (e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
}

/*
  Use neighbor list
  Use a table for Energy
*/
void single_flip_table_list(double beta) {
  for (int i = 0; i < N; i++) {
    int ns = 0;
    ns += spin[neighbor[i][0]];
    ns += spin[neighbor[i][1]];
    ns += spin[neighbor[i][2]];
    ns += spin[neighbor[i][3]];
    ns *= spin[i];
    if (ns < 0 || e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
}

/*
 Calculate bulk without adjusting PBC
*/

void single_flip_bulk(double beta) {

  // Bulk
  for (int iy = 1; iy < L - 1; iy++) {
    for (int ix = 1; ix < L - 1; ix++) {
      int i = ix + iy * L;
      int ns = 0;
      ns += spin[i + 1];
      ns += spin[i - 1];
      ns += spin[i + L];
      ns += spin[i - L];
      ns *= spin[i];
      if (ns < 0 || e_table[ns / 2 + 2] > ud(mt)) {
        spin[i] = -spin[i];
      }
    }
  }

  // site_i = (j, 0)
  for (int j = 0; j < L; j++) {
    int i = j;
    int ns = 0;
    ns += spin[neighbor[i][0]];
    ns += spin[neighbor[i][1]];
    ns += spin[neighbor[i][2]];
    ns += spin[neighbor[i][3]];
    ns *= spin[i];
    if (ns < 0 || e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }

  // site_i = (j, L-1)
  for (int j = 0; j < L; j++) {
    int i = j + (L - 1) * L;
    int ns = 0;
    ns += spin[neighbor[i][0]];
    ns += spin[neighbor[i][1]];
    ns += spin[neighbor[i][2]];
    ns += spin[neighbor[i][3]];
    ns *= spin[i];
    if (ns < 0 || e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
  // site_i = (0, j)
  for (int j = 1; j < L - 1; j++) {
    int i = j * L;
    int ns = 0;
    ns += spin[neighbor[i][0]];
    ns += spin[neighbor[i][1]];
    ns += spin[neighbor[i][2]];
    ns += spin[neighbor[i][3]];
    ns *= spin[i];
    if (ns < 0 || e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }

  // site_i = (L-1, j)
  for (int j = 1; j < L - 1; j++) {
    int i = L - 1 + j * L;
    int ns = 0;
    ns += spin[neighbor[i][0]];
    ns += spin[neighbor[i][1]];
    ns += spin[neighbor[i][2]];
    ns += spin[neighbor[i][3]];
    ns *= spin[i];
    if (ns < 0 || e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
}

void single_flip(double beta) {
  //single_flip_table_list(beta);
  //single_flip_table_mod_tableonly(beta);
  single_flip_bulk(beta);
}

/*
Order Parameter (m^2)
*/
double magnetization_square(void) {
  double m = std::accumulate(spin.begin(), spin.end(), 0.0);
  m /= static_cast<double>(N);
  return m * m;
}

void make_table(double beta) {
  for (int i = 0; i < 5; i++) {
    int ns = (i - 2) * 4;
    e_table[i] = exp(-ns * beta);
  }
}

void domc(void) {
  const int T_LOOP = 10000; // Thermalization loop
  const int O_LOOP = 1000;  // Number of observations
  const int ND = 20;        // Number of temperatures to be observed
  double bs = 0.3;          //Start inverse temperature  (The highest temperature)
  double be = 0.6;          //End inverse temperature (The lowest temperature)
  for (int i = 0; i < ND; i++) {
    double beta = bs + (be - bs) * (double)i / (double)ND;
    make_table(beta);
    std::cerr << "beta = " << beta << std::endl;
    std::fill(spin.begin(), spin.end(), 1);
    for (int j = 0; j < T_LOOP; j++) {
      single_flip(beta);
    }
    double sm2 = 0.0;
    for (int k = 0; k < O_LOOP; k++) {
      single_flip(beta);
      double m2 = magnetization_square();
      sm2 += m2;
    }
    sm2 /= static_cast<double>(O_LOOP);

    std::cout << beta << " " << sm2 << std::endl;
  }
}

/*
  Check the update speed at the critical point beta_c = log(1+sqrt(2))/2
 */
void speed_test(const char *name, void (*single_flip)(double)) {
  double beta_c = 0.5 * log(1.0 + sqrt(2.0));
  make_table(beta_c);
  std::fill(spin.begin(), spin.end(), 1);
  mt.seed(1);
  const int LOOP = 100000;
  const auto s = std::chrono::system_clock::now();
  for (int i = 0; i < LOOP; i++) {
    single_flip(beta_c);
  }
  const auto e = std::chrono::system_clock::now();
  const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();
  std::cout << std::setprecision(4);
  std::cout << name << "\t" << magnetization_square() << "\t" << elapsed << "[msec]" << std::endl;
}

int main(void) {
  init_neighbors();
#ifdef MC
  domc();
#else
  speed_test("calc      + if  ", single_flip_calc_if);
  speed_test("calc      + mod ", single_flip_calc_mod);
  speed_test("table     + mod ", single_flip_table_mod);
  speed_test("tableonly + mod ", single_flip_table_mod_tableonly);
  speed_test("table     + list", single_flip_table_list);
  speed_test("bulk + border   ", single_flip_bulk);
#endif
}
