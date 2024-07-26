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
  single_flip_table_list(beta);
  //single_flip_table_mod_tableonly(beta);
  //single_flip_bulk(beta);
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

void relaxation(double beta) {
  make_table(beta);
  std::cerr << "beta = " << beta << std::endl;
  double sm2 = 0.0;
  const int O_LOOP = 1000;
  const int N_SAMPLE = 1000;
  std::vector<double> data(1000, 0.0);
  for (int j = 0; j < N_SAMPLE; j++) {
    std::fill(spin.begin(), spin.end(), 1);
    for (int i = 0; i < O_LOOP; i++) {
      single_flip(beta);
      double m2 = magnetization_square();
      data[i] += m2;
    }
  }
  for (int i = 0; i < O_LOOP; i++) {
    std::cout << i << " " << data[i] / N_SAMPLE << std::endl;
  }
}

int main(void) {
  init_neighbors();
  relaxation(0.42);
}
