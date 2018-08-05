Autocorrelation Function of Magnetization of 2D Ising Model
===

# Summary

Sample source code of Monte Calro simulation for the two-dimensional
Ising model on the square lattice.
Autocorrelation functions of magnetization at the critical temperature
are calculated by three methods, Single-Flip, Swendsen-Wang, and
Wolff algorithm, respectively.

# Usage

    $ make clean
    $ make acf


# Files

- main.cc

  A source file
- single.dat , sw.dat, wolff.dat

  Time evolutions of square of magnetization at the critical point.

- single.acf, sw.acf, wolff.acf

  Autocorrelation function of the magetization.

- auto.rb

  Ruby scripts to calculate autocorrelation function.

- acf.plt

  Plot file for gnuplot

- acf.png

  Plot image by gnuplot
