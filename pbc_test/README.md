# Periodic Boundary Adjustment

## Summary

I investigated the efficient implementation of the periodic boundary condition.

## Benchmark results

### `calc + if` (Baseline)

Consider the ferromagnetic Ising model on the LxL square lattice. We adopt the Metropolis single-spin udpate algorithm. The program would be something like as follows.

```cpp
  for (int i = 0; i < N; i++) {
    int ns = 0;
    ns += spin[i + 1];
    ns += spin[i - 1];
    ns += spin[i + L];
    ns += spin[i - L];
    ns *= spin[i];
    if (ns < 0 || exp(-2.0 * ns * beta) > ud(mt)) {
      spin[i] = -spin[i];
    }
  }
```

Since the above program does not care the boundary condition, it touches the memory out of boundary.

If we adopt the perodic boundary condition, we have to adjust the boundary like as,

```cpp
  if (ix < 0) ix += L;
  if (ix >= L) ix -= L;
  if (iy < 0) iy += L;
  if (iy >= L) iy -= L;
```

So, the program would be,

```cpp
int pos2index_if(int ix, int iy) {
  if (ix < 0) ix += L;
  if (ix >= L) ix -= L;
  if (iy < 0) iy += L;
  if (iy >= L) iy -= L;
  return ix + iy * L;
}

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
```

We refer the above program to `calc+if` and adopt it as the baseline. We consider the square system with L=64. We updates 100000 times at the critical point. At some machine, it took 16151 [msec].

### `calc + mod`

One can find that the function `pos2index_if` is inefficient. We can use `mod` instead of `if` as,

```cpp
int pos2index_mod(int ix, int iy) {
  ix = (ix + L) % L;
  iy = (iy + L) % L;
  return ix + iy * L;
}

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
```

The elapse time becomes 16064[msec]. The efficiency was not improved so much.

### `table + mod`

While we calculate `exp` funciton every time, the number of possible energy is limited. So we adopt a table for the energy difference.

```cpp
void make_table(double beta) {
  for (int i = 0; i < 5; i++) {
    int ns = (i - 2) * 4;
    e_table[i] = exp(-ns * beta);
  }
}

```cpp
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
```

We have 7488 [msec], we achieves 2x speedup.

### `tableonly + mod`

When `ns < 0`, then `exp(-2.0 * ns * beta)` is always larger than 1. Therefore, we can replace

```cpp
    if (ns < 0 || exp(-2.0 * ns * beta) > ud(mt)) {
```

with

```cpp
    if (e_table[ns / 2 + 2] > ud(mt)) {
```

The result was 8009 [msec], so the performance was degraded.

### `table + list`

We construct the neighborlist of spins like as,

```cpp
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
```

Using the neighborlist, we can implement spin updates as

```cpp
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
```

The result is 6815 [msec], so performance was slightly improved.

### `bulk + border`

For the spins in the bulk, i.e., 1 < x < L-1 and 1 < y < L-1, we do not have to take care of the boundary conditions.

So we can updates spins in the bulk and spins at the border separately. The code is,

```cpp
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
```

The result was 6595 [msec], which is the fastest.

## Results

`name`, `magnetization`, `elapsed time` are listed.

### Linux + GCC

* g++ (GCC) 7.2.0
* Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz
* `-O3 -std=c++11 -march=native`
* Red Hat Enterprise Linux Server release 7.4 (Maipo)

```txt
calc      + if  	0.3537	16151[msec]
calc      + mod 	0.3537	16064[msec]
table     + mod 	0.3537	7488[msec]
tableonly + mod 	0.2853	8009[msec]
table     + list	0.3537	6815[msec]
bulk + border   	0.5839	6595[msec]
```

### Linux + Intel Compiler

* icpc (ICC) 18.0.5 20180823
* Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz
* `-O3 -std=c++11 -march=native`
* Red Hat Enterprise Linux Server release 7.4 (Maipo)

```txt
calc      + if  	0.3537	33342[msec]
calc      + mod 	0.3537	32431[msec]
table     + mod 	0.3537	6152[msec]
tableonly + mod 	0.2853	6458[msec]
table     + list	0.3537	30158[msec]
bulk + border   	0.5839	29927[msec]
```

### Mac + Clang

* Apple clang version 11.0.0 (clang-1100.0.33.8)
* Intel Core i5 3.3 GHz
* macOS Mojave 10.14.6
* `-O3 -std=c++11 -march=native`

```txt
calc      + if  	0.3537	8778[msec]
calc      + mod 	0.3537	7553[msec]
table     + mod 	0.3537	5622[msec]
tableonly + mod 	0.2853	5980[msec]
table     + list	0.3537	5574[msec]
bulk + border   	0.5839	5398[msec]
```