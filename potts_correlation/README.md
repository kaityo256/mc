------------------------------------------------------------------------
# Calculate pair correlation function of Q-state Potts model
------------------------------------------------------------------------

## Summary

Sample code to calculate pair correlation function of Q-state Potts model on the square lattice. You need FFTW3 library to build this.

## Usage

    $ make
    $ ./potts2d > test.dat

## Parameters

- Q: Number of states (when Q=2, it is called Ising model)
- LX, LY: System Size
- T\_LOOP: A number of MC steps for thermalization
- O\_LOOP: A number of MC steps for observation
