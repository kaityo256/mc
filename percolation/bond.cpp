// Bond Percolation

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

int L; // System Size
int N; // Number of Sites

std::vector<int> parent;

void init(int size) {
  L = size;
  N = L * L;
  parent.resize(N);
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
  if (ud(mt) > p) return;
  unite(i, j);
}

int pos2index(int ix, int iy) {
  ix = (ix + L) % L;
  iy = (iy + L) % L;
  return ix + iy * L;
}

void one_step(double p, std::mt19937 &mt) {
  for (int i = 0; i < N; i++) {
    parent[i] = i;
  }
  for (int iy = 0; iy < L - 1; iy++) {
    for (int ix = 0; ix < L - 1; ix++) {
      int i = pos2index(ix, iy);
      connect(i, pos2index(ix + 1, iy), p, mt);
      connect(i, pos2index(ix, iy + 1), p, mt);
    }
  }
}

double crossing_probability(void) {
  for (int ix1 = 0; ix1 < L; ix1++) {
    int i = find(pos2index(ix1, 0));
    int ci = find(i);
    for (int ix2 = 0; ix2 < L; ix2++) {
      int j = find(pos2index(ix2, L - 1));
      int cj = find(j);
      if (ci == cj)
        return 1.0;
    }
  }
  return 0.0;
}

double percolation_probability(void) {
  std::vector<int> size(N, 0);
  for (int i = 0; i < N; i++) {
    int ci = find(i);
    size[ci]++;
  }
  int max = *std::max_element(size.begin(), size.end());
  return static_cast<double>(max) / N;
}

void mc(double p) {
  const int observe_loop = 1000;
  std::mt19937 mt;
  double cp = 0;
  double pp = 0.0;
  for (int i = 0; i < observe_loop; i++) {
    one_step(p, mt);
    cp += crossing_probability();
    pp += percolation_probability();
  }
  cp /= static_cast<double>(observe_loop);
  pp /= static_cast<double>(observe_loop);
  std::cout << p << " " << cp << " " << pp << std::endl;
}

int main(void) {
  int ND = 50;
  int size = 32;
  init(size);
  for (int i = 0; i <= ND; i++) {
    double p = static_cast<double>(i) / ND;
    mc(p);
  }
}
