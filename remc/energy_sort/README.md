------------------------------------------------------------------------
#Adjustment of a set of temperatures for Replica Exchange MC
------------------------------------------------------------------------
##Summary

Sample source code of Replica Exchange Monte Carlo (REMC).
Each replica has fixed energy. After several steps of REMC,
the replicas are sorted by their energy, i.e, the replica
with higher energy has higher temperature.

-----------------------------------------------------------------------
##Usage

   $ make clean
   $ make

-----------------------------------------------------------------------
##Files

-exmc.rb
  A source file

-index.dat  
  Indeces of replicas. 

-index.plt
  Plot file for gnuplot

-index.png
  Plot image by gnuplot

------------------------------------------------------------------------
