# PDA_C
A prototype of a serial object-oriented C++ implementation of the Population Dynamics Algorithm
that was published by Charlebois et al. in 2011 in Communications in Computational Physics.

Notes: This code has not been optimized or tested extensively against analytic solutions. Segmentation fault occurs after ~20k cells. This algorithm must be compiled with "mtrand.cpp" - a C++ implementation of the 
Mersenne twister pseudorandom number generator.

*For a parallel C++ implementation of the PDA see the "gepda" repository of GitHub user "yuchenhou". 
