# 2-Dimensional Finite Difference Hartree-Fock program for diatomic molecules
                                                                            
This program finds virtually exact solutions of the Hartree-Fock and
density functional theory equations for diatomic molecules and atoms
(the quality of a solution depends on grid size and arithmetic
precision used). The lowest energy eigenstates of a given irreducible
representation and spin can be obtained.

The program can be also be used to obtain the ground and excited
states of one-electron systems with the (smoothed) Coulomb and
Kramers-Hennenberger potentials.

Single particle two-dimensional numerical functions (orbitals) are
used to construct an antisymmetric many-electron wave function of the
restricted open-shell Hartree-Fock model. The orbitals are obtained by
solving the Hartree-Fock equations in the form of the coupled
two-dimensional second-order (elliptic) partial differential equations
(PDE). The Coulomb and exchange potentials are obtained as solutions
of the corresponding Poisson equations. The PDEs are disretized by the
8th-order central difference stencil on a two-dimensional grid and the
resulting large and sparse system of linear equations is solved by the
(multicolour) successive overrelaxation method ((MC)SOR). The
self-consistent-field iterations are interwoven with the (MC)SOR ones
and orbital energies and normalization factors are used to monitor the
convergence. The accuracy of solutions depends mainly on the grid and
the system under consideration.

See the following articles for the detailed description of the program
and examples of its usage and accuracy:

* L. Laaksonen, P. Pyykkö, and D. Sundholm, Fully numerical Hartree-Fock methods for molecules, Comp. Phys. Reports 4 (1986) 313-344. http://doi.org/10.1016/0167-7977(86)90021-3
* L. Laaksonen, D. Sundholm, and P. Pyykkö, in "Scientific Computing in Finland", Eds. K. Kankaala and R. Nieminen, Research Report R1/89, Centre for Scientific Computing, (1989) p. 183.
* P. Pyykkö, in Numerical Determination of the Electronic Structure of Atoms, Diatomic and Polyatomic Molecules (NATO ASI Series C271) Eds. M. Defranceschi and J. Delhalle, (1989) p. 161.     http://doi.org/10.1007/978-94-009-2329-4
* J. Kobus, Finite-difference versus finite-element methods, Chem. Phys. Lett. 202 (1993) 7-12. http://doi.org/10.1016/0009-2614(93)85342-L
* J. Kobus, Vectorizable algorithm for the (multicolour) successive overrelaxation method, Comput. Phys. Commun. 78 (1994) 247-255. http://doi.org/10.1016/0010-4655(94)90003-5
* J. Kobus, L. Laaksonen, D. Sundholm, A numerical Hartree-Fock program for diatomic molecules, Comp. Phys. Commun. 98 (1996) 346-358. http://doi.org/10.1016/0010-4655(96)00098-7
* J. Kobus, Numerical Hartree-Fock methods for diatomic molecules, Handbook of Molecular Physics and Quantum Chemistry (Chichester), ed. S. Wilson (Wiley, 2002)
* J. Kobus, Hartree-Fock limit values of multipole moments, polarizabilities and hyperpolarizabilities for atoms and diatomic molecules, Comp. Lett. 3 (2007) 71-113. http://doi.org/10.1163/157404007782913408
* J. Kobus, A finite difference Hartree-Fock program for atoms and diatomic molecules, Comp. Phys. Commun. 184 (2013) 799-811. http://dx.doi.org/10.1016/j.cpc.2012.09.033

The last paper in the above list has an up-to-date description of the
operating principles and layout of the program.

The programming language used is Fortran 90. Only a Fortran compiler
is necessary to compile the program. The program can be built using
the either CMake, or autotools.

To build x2dhf with CMake, first configure the build with `cmake
.`, after which you can build the program with `$ make`.

To build x2dhf with autotools, run `$ ./autoconf.sh` to generate the
configure script, after which run configure `$ ./configure` and make
the program with `$ make`.

The file doc/users-guide.pdf contains the description of the program's
input data and examples of its usage. Examples are included in
examples/

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.
                                                                      
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
