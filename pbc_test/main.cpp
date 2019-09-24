#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

const int L = 64;    // System Size
const int N = L * L; // Number of spins
double e_table[5];

std::vector<int> spin(N), cluster(N), flip(N);
std::vector<std::array<int, 4>> f_neighbor(N);

std::mt19937 mt(1);
std::uniform_real_distribution<double> ud(0.0, 1.0);

class svariable {
private:
  std::vector<double> data, data2;

public:
  void add(double v) {
    data.push_back(v);
    data2.push_back(v * v);
  }

  double average() {
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
  }

  double error() {
    double ave = average();
    double ave2 = std::accumulate(data2.begin(), data2.end(), 0.0) / data2.size();
    double v = (ave2 - ave * ave) / (data.size() - 1);
    return sqrt(v);
  }
};

int pos2index(int ix, int iy) {
  ix = (ix + L) % L;
  iy = (iy + L) % L;
  return ix + iy * L;
}

void single_flip(double beta) {
  for (int i = 0; i < N; i++) {
    int ns = 0;
    ns += spin[f_neighbor[i][0]];
    ns += spin[f_neighbor[i][1]];
    ns += spin[f_neighbor[i][2]];
    ns += spin[f_neighbor[i][3]];
    ns *= spin[i];
    if (e_table[ns / 2 + 2] > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
}

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
  const int T_LOOP = 1000; // Thermalization loop
  const int O_LOOP = 1000; // Number of observations
  const int ND = 20;       // Number of temperatures to be observed
  double bs = 0.3;         //Start inverse temperature  (The highest temperature)
  double be = 0.6;         //End inverse temperature (The lowest temperature)
  for (int i = 0; i < ND; i++) {
    double beta = bs + (be - bs) * (double)i / (double)ND;
    make_table(beta);
    std::cerr << "beta = " << beta << std::endl;
    std::fill(spin.begin(), spin.end(), 1);
    for (int j = 0; j < T_LOOP; j++) {
      single_flip(beta);
    }
    svariable sm2;
    svariable sm4;
    for (int k = 0; k < O_LOOP; k++) {
      single_flip(beta);
      double m2 = magnetization_square();
      sm2.add(m2);
      sm4.add(m2 * m2);
    }
    double am2 = sm2.average();
    double am4 = sm4.average();
    double U = am4 / (am2 * am2);

    std::cout << beta << " " << sm2.average() << " " << sm2.error();
    std::cout << " " << U << std::endl;
  }
}

void init_neighbors() {
  for (int iy = 0; iy < L; iy++) {
    for (int ix = 0; ix < L; ix++) {
      int i = pos2index(ix, iy);
      f_neighbor[i][0] = pos2index(ix + 1, iy);
      f_neighbor[i][1] = pos2index(ix - 1, iy);
      f_neighbor[i][2] = pos2index(ix, iy + 1);
      f_neighbor[i][3] = pos2index(ix, iy - 1);
    }
  }
}

/*
  Check the update speed at the critical point beta_c = log(1+sqrt(2))/2
 */
void speed_test(void) {
  double beta_c = 0.5 * log(1.0 + sqrt(2.0));
  make_table(beta_c);
  std::fill(spin.begin(), spin.end(), 1);
  const int LOOP = 10000;
  const auto s = std::chrono::system_clock::now();
  for (int i = 0; i < LOOP; i++) {
    single_flip(0.001);
  }
  const auto e = std::chrono::system_clock::now();
  const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();
  std::cout << beta_c << " " << magnetization_square() << " " << elapsed << "[msec]" << std::endl;
}

int main(void) {
  init_neighbors();
  //domc();
  speed_test();
}
