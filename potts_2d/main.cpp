#include <algorithm>
#include <array>
#include <iostream>
#include <random>
#include <vector>

int L; // System Size
int N; // Number of Spins
int Q; // Number of States

const int thermalize_loop = 100;
const int observe_loop = 10000;

std::vector<int> spins, newspins;
std::vector<std::array<int, 2>> neighbor;
std::vector<int> parent;

int pos2index(int ix, int iy) {
  ix = (ix + L) % L;
  iy = (iy + L) % L;
  return ix + iy * L;
}

void init(int size) {
  L = size;
  N = L * L;
  spins.resize(N);
  newspins.resize(N);
  neighbor.resize(N);
  parent.resize(N);
  // Initialize Neighbors
  for (int iy = 0; iy < L; iy++) {
    for (int ix = 0; ix < L; ix++) {
      int i = pos2index(ix, iy);
      neighbor[i][0] = pos2index(ix + 1, iy);
      neighbor[i][1] = pos2index(ix, iy + 1);
    }
  }
}

int find(int i) {
  while (i != parent[i]) {
    i = parent[i];
  }
  return i;
}

void unite(int i, int j) {
  i = find(i);
  j = find(j);
  parent[j] = i;
}

void connect(int i, int j, double p, std::mt19937 &mt) {
  std::uniform_real_distribution<> ud(0.0, 1.0);
  if (spins[i] != spins[j]) return;
  if (ud(mt) > p) return;
  unite(i, j);
}

void one_step(double beta, std::mt19937 &mt) {
  std::uniform_int_distribution<> ud(0, Q - 1);
  for (int i = 0; i < N; i++) {
    parent[i] = i;
  }
  // Connect
  const double p = 1.0 - exp(-beta);
  for (int i = 0; i < N; i++) {
    connect(i, neighbor[i][0], p, mt);
    connect(i, neighbor[i][1], p, mt);
  }
  // Flip
  for (int i = 0; i < N; i++) {
    newspins[i] = ud(mt);
  }
  for (int i = 0; i < N; i++) {
    spins[i] = newspins[find(i)];
  }
}

double magnetization(void) {
  std::vector<double> m(Q, 0.0);
  for (int i = 0; i < N; i++) {
    m[spins[i]] += 1;
  }
  for (int i = 0; i < Q; i++) {
    m[i] /= static_cast<double>(N);
  }
  double m2 = 0.0;
  for (int i = 0; i < Q; i++) {
    m2 += m[i] * m[i];
  }
  for (int i = 0; i < Q - 1; i++) {
    for (int j = i + 1; j < Q; j++) {
      m2 -= 2.0 * m[i] * m[j] / (Q - 1);
    }
  }
  return m2;
}

void mc(double beta) {
  std::mt19937 mt;
  for (int i = 0; i < thermalize_loop; i++) {
    one_step(beta, mt);
  }
  double sm2 = 0.0;
  double sm4 = 0.0;
  for (int i = 0; i < observe_loop; i++) {
    one_step(beta, mt);
    double m2 = magnetization();
    sm2 += m2;
    sm4 += m2 * m2;
  }
  sm2 /= static_cast<double>(observe_loop);
  sm4 /= static_cast<double>(observe_loop);
  double u = sm4 / sm2 / sm2;
  std::cout << beta << " " << sm2 << " " << u << std::endl;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "./a.out Q systemsize beta_start beta_end" << std::endl;
    return 0;
  }
  Q = std::stoi(argv[1]);
  L = std::stoi(argv[2]);
  init(L);
  double bs = std::stof(argv[3]);
  double be = std::stof(argv[4]);
  const int ND = 40;
  for (int i = 0; i < ND; i++) {
    double beta = bs + (be - bs) * i / ND;
    mc(beta);
  }
}