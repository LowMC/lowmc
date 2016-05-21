# LowMC implementation
This is a C++ implementation of the LowMC block cipher family. The
parameters (block size, number of S-boxes in the substitution layer, number
of rounds, key size) are defined in `LowMC.h`. Compilation requires support
of C++11 features.

The files `LowMC.h` and `LowMC.cpp` contain the relevant code. `test.cpp`
contains a short example usage.
