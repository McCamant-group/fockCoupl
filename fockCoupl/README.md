This program takes NxN Fock and overlap matrices in the atomic orbital basis, where N is the total number of basis functions, and diaganalizes user defined blocks which refer to donor and acceptor partitions of the system.

The program is currently written so that a is the acceptor and b is the donor. In the future this will be generalized so that the order does not matter. The input Fock and Overlap matrices need to come from an input structure ordered such that the first set of basis functions will all correspond to the acceptor (a) and the second set of basis functions will all correspond to the donor (b). Future update will include exctraction script to parse either ORCA or Gaussian output files for the Fock and overlap matrices.

The user must input the number of basis functions associated with partition a and partition b.

Once the Fock matrix is block diagonalized, the eigenvectors are used to calculate the off-diagonal electronic coupling matrix elements between a and b. 

This coupling is then used to compute electron transfer rates from a selected orbital on site b to all orbitals on site a. This calculation utilizes the electronic coupling from the donor orbital to all acceptor orbitals as well as their site energy differences. The reorganization energy needs to be input by the user. 

The program will output plots of the most relevant data including the b_orbital coupling to every a acceptor orbital, the partitioned eigenvalues, the eT rate from b_orbital to every a acceptor orbital, and the time constant and population dynamics which result from the calculation. These results are also printed to output files for analysis in other programs. 