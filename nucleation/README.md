Simple Nucleation 
===

# Summary

Sample source code of Monte Carlo simulation for a random walk in
a potential, which simulates the nucleation of a droplet in
supersaturated vapor. The potential function is 

$$
  V(x) = x^2 - x^3,
$$

where $x$ is a radius of a droplet.
First passage time (nucleation time) is calculated as a function
of inverse temperature.


# Usage

    $ make clean
    $ make


# Files

- nucleation.rb

  A source file

- nucleation.dat

  Nucleation Time vs. Inverse temperature

- nucleation.plt

  Plot file for gnuplot

- nucleation.png

  Plot image by gnuplot
