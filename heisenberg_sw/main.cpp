#include <iostream>
#include <cstdlib>
#include <cmath>
#include <random>
#include <algorithm>

const int L = 32;         // System Size
const int T_LOOP = 10000; // Thermalization loop
const int O_LOOP = 5000;  // Number of average
const int ND = 20;        // Number of temperatures to be observed

class Spin {
  public:
    double sx, sy, sz;
    Spin() {
      sx = 0.0;
      sy = 0.0;
      sz = 1.0;
    };
    Spin(double x, double y, double z) {
      sx = x;
      sy = y;
      sz = z;
    };

    Spin(double z, double phi) {
      const double c = sqrt(1.0 - z * z);
      sx = c * cos(phi);
      sy = c * sin(phi);
      sz = z;
    };

    Spin operator+(const Spin &s) {
      Spin s2;
      s2.sx = sx + s.sx;
      s2.sy = sy + s.sy;
      s2.sz = sz + s.sz;
      return s2;
    };

    Spin operator-(const Spin &s) {
      Spin s2;
      s2.sx = sx - s.sx;
      s2.sy = sy - s.sy;
      s2.sz = sz - s.sz;
      return s2;
    };

    void operator+=(const Spin &s) {
      sx += s.sx;
      sy += s.sy;
      sz += s.sz;
    };

    void operator-=(const Spin &s) {
      sx -= s.sx;
      sy -= s.sy;
      sz -= s.sz;
    };

    Spin operator*(const double c) const {
      Spin s2(sx, sy, sz);
      s2.sx *= c;
      s2.sy *= c;
      s2.sz *= c;
      return s2;
    }

    double operator*(const Spin &s) const {
      return sx * s.sx + sy * s.sy + sz * s.sz;
    }

    void operator=(const Spin &s) {
      sx = s.sx;
      sy = s.sy;
      sz = s.sz;
    };

    void operator/=(const double d) {
      sx /= d;
      sy /= d;
      sz /= d;
    };

    Spin flip(const Spin &r) const {
      Spin s(sx, sy, sz);
      return s - (r * (r * 2.0 * (*this)));
    };

    const double square(void) const {
      return sx * sx + sy * sy + sz * sz;
    };

};

Spin spins[L][L];
int cluster[L * L];
bool flip[L * L];

std::mt19937 mt(1);
std::uniform_real_distribution<double> ud(0.0, 1.0);

void init(void) {
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      spins[i][j].sx = 0.0;
      spins[i][j].sy = 0.0;
      spins[i][j].sz = 1.0;
    }
  }
}

void single_flip(const double beta) {
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      Spin ms(0.0, 0.0, 0.0);
      if (ix != 0) ms += spins[ix - 1][iy];
      else ms += spins[L - 1][iy];
      if (iy != 0) ms += spins[ix][iy - 1];
      else ms += spins[ix][L - 1];
      if (ix != L - 1) ms += spins[ix + 1][iy];
      else ms += spins[0][iy];
      if (iy != L - 1) ms += spins[ix][iy + 1];
      else ms += spins[ix][0];
      Spin ns(ud(mt) * 2.0 - 1, ud(mt) * 2.0 * M_PI);
      Spin ds = spins[ix][iy] - ns;
      const double de = ds * ms;
      if (de < 0 || exp(-beta * de) > ud(mt)) {
        spins[ix][iy] = ns;
      }
    }
  }
}

int get_cluster_number(int index) {
  int i = index;
  while (i != cluster[i]) {
    i = cluster[i];
  }
  return i;
}

void connect_bond(int ix1, int iy1, int ix2, int iy2) {
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

void cluster_flip(double beta) {
  for (int i = 0; i < L * L; i++) {
    cluster[i] = i;
    flip[i] = (ud(mt) > 0.5);
  }
  const Spin r(ud(mt) * 2.0 - 1, ud(mt) * 2.0 * M_PI);

  for (int ix = 0; ix < L - 1; ix++) {
    for (int iy = 0; iy < L; iy++) {
      const double dE = 2.0 * (r * spins[ix][iy]) * (r * spins[ix + 1][iy]);
      if (ud(mt) < 1.0 - exp(-beta * dE)) {
        connect_bond(ix, iy, ix + 1, iy);
      }
    }
  }

  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L - 1; iy++) {
      const double dE = 2.0 * (r * spins[ix][iy]) * (r * spins[ix][iy + 1]);
      if (ud(mt) < 1.0 - exp(-beta * dE)) {
        connect_bond(ix, iy, ix, iy + 1);
      }
    }
  }

  for (int iy = 0; iy < L; iy++) {
    const double dE = 2.0 * (r * spins[L - 1][iy]) * (r * spins[0][iy]);
    if (ud(mt) < 1.0 - exp(-beta * dE)) {
      connect_bond(L - 1, iy, 0, iy);
    }
  }

  for (int ix = 0; ix < L; ix++) {
    const double dE = 2.0 * (r * spins[ix][L - 1]) * (r * spins[ix][0]);
    if (ud(mt) < 1.0 - exp(-beta * dE)) {
      connect_bond(ix, L - 1, ix, 0);
    }
  }
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      const int c = get_cluster_number(ix + iy * L);
      if (flip[c]) {
        spins[ix][iy] = spins[ix][iy].flip(r);
      }
    }
  }
}

double magnetization_square(void) {
  Spin s(0.0, 0.0, 0.0);
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      s += spins[ix][iy];
    }
  }
  s /= static_cast<double>(L * L);
  return s.square();
}

double energy(void) {
  double e = 0.0;
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      Spin s = spins[ix][iy];
      if (ix != 0) e -= s * spins[ix - 1][iy];
      else e -= s * spins[L - 1][iy];
      if (iy != 0) e -= s * spins[ix][iy - 1];
      else e -= s * spins[ix][L - 1];
    }
  }
  return e /= static_cast<double>(L * L);
}

void mc_onestep(double beta) {
#ifdef CLUSTER
  cluster_flip(beta);
#else
  single_flip(beta);
#endif
}

void domc(void) {
  double ts = 0.01;
  double te = 3.0;
  for (int i = 0; i < ND; i++) {
    double beta = 1.0 / (ts + (te - ts) * (double)i / (double)ND);
    fprintf(stderr, "T = %f\n", 1.0 / beta);
    init();
    for (int j = 0; j < T_LOOP; j++) {
      mc_onestep(beta);
    }
    double m2_sum = 0.0;
    double m4_sum = 0.0;
    double e_sum = 0.0;
    double e2_sum = 0.0;
    for (int j = 0; j < O_LOOP; j++) {
      mc_onestep(beta);
      const double m2 = magnetization_square();
      m2_sum += m2;
      m4_sum += m2 * m2;
      const double e = energy();
      e_sum += e;
      e2_sum += (e * e);
    }
    m2_sum /= O_LOOP;
    m4_sum /= O_LOOP;
    e_sum /= O_LOOP;
    e2_sum /= O_LOOP;
    const double c = L * L * (e2_sum -  e_sum * e_sum) * beta * beta;
    const double binder = m4_sum / (m2_sum * m2_sum);
    printf("%f %f %f %f\n", 1.0 / beta, m2_sum, c, binder);
  }
}

int main(void) {
#ifdef CLUSTER
  fprintf(stderr, "Swendsen-Wang\n");
  printf("# beta, Magnetization^2, Specific Heat, Binder Parameter\n");
#else
  fprintf(stderr, "Single Flip\n");
  printf("# beta, Magnetization^2, Specific Heat, Binder Parameter\n");
#endif
  domc();
}
