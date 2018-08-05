Adjustment of a set of temperatures for Replica Exchange MC
===

# Summary

Sample source code of Replica Exchange Monte Carlo (REMC) with 
adjustment of a set of temperatures.
If the specific heat of the system is constant, i.e.,
E(beta) \sim 1/beta, then the temperatures of the replicas
should be geometric progression in order to achieve
uniform exchange rates between adjacent temperatures.

In this sample, the exchante rate is adjusted to be 0.6.

# Usage

    $ make clean
    $ make

# Files

* [exmc.rb](exmc.rb)    
  A source file

* [beta.dat](beta.dat)    
  Tempeartures obtained by REMC

* [beta.plt](beta.plt)    
  Plot file for gnuplot

* [beta.png](beta.png)    
  Plot image by gnuplot
