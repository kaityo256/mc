------------------------------------------------------------------------
#Adjustment of a set of temperatures for Replica Exchange MC
------------------------------------------------------------------------
##Summary

Sample source code of Replica Exchange Monte Carlo (REMC) with 
adjustment of a set of temperatures.
If the specific heat of the system is constant, i.e.,
E(\beta) \sim 1/\beta, then the temperatures of the replicas
should be geometric progression in order to achieve
uniform exchange rates between adjacent temperatures.

-----------------------------------------------------------------------
#Usage

    $ make clean
    $ make

-----------------------------------------------------------------------
#Files

-exmc.rb
  A source file

-beta.dat  
  Tempeartures obtained by REMC

-beta.plt
  Plot file for gnuplot

-beta.png
  Plot image by gnuplot

------------------------------------------------------------------------