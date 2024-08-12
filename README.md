# Hartree-Fock---Square-lattice
This project aims to solve the Hubbard model for a square lattice using the Hartree-Fock approximation. 
By varying the number of electrons in the lattice and the Coulombian repulsion, we have constructed the magnetic phase diagram, 
which is in agreement with the literature (Y. Claveau et al., 2014, Eur. J. Phys., 35, 035023).

The basic code was developed using Fortran90, but the presentation of the results used Python.

Main code:
 - Makefile
 - interface.f90
 - run2.f90 (Converges the magnetic phases)

Results code:
 - ler_phase.py

Results files:
 - ANTI_ne1.txt ; ANTI_ne2.txt ; ANTI_U1.txt ; ANTI_U2.txt (Antiferromagnetic phase)
 - FERRO_ne1.txt ; FERRO_ne2.txt ; FERRO_U1.txt ; FERRO_U2.txt (Ferromagnetic phase)
 - PARA_ne1.txt ; PARA_ne2.txt ; PARA_U1.txt ; PARA_U2.txt (Paramagnetic phase)

Results figure:
 - phase_diagram.pdf
