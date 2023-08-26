Two-Dimensional Finite Difference Hartree-Fock program for diatomic molecules (version 3.0)
-------------------------------------------------------------------------------------------------  

                                                                            
This program finds virtually exact solutions of the Hartree-Fock and density
functional theory equations for diatomic molecules and atoms (the quality of a
solution depends on grid size and arithmetic precision used). The lowest energy
eigenstates of a given irreducible representation and spin can be obtained. The DFT
support is facilitated via the libxc library (see `-l|L` options of the `./x2dhfctl`
script). `lxcctl` script is provided to help test the support of the `x2dhf` program
for the libxc functionals.

The program can also used to obtain the ground and excited states of one-electron
systems with the (smoothed) Coulomb and Kramers-Hennenberger potentials.

Single particle two-dimensional numerical functions (orbitals) are used to construct
an anti-symmetric many-electron wave function of the restricted open-shell
Hartree-Fock model. The orbitals are obtained by solving the Hartree-Fock equations
in the form of the coupled two-dimensional second-order (elliptic) partial
differential equations (PDE). The Coulomb and exchange potentials are obtained as
solutions of the corresponding Poisson equations. The PDEs are discretized by the
8th-order central difference stencil on a two-dimensional grid. The resulting large
and sparse system of linear equations is solved by the (multi-colour) successive
overrelaxation method ((MC)SOR). The self-consistent-field iterations are interwoven
with the (MC)SOR ones and orbital energies and normalisation factors are used to
monitor the convergence. The accuracy of solutions depends mainly on the grid-size and the
system under consideration.

See the following articles for the detailed description of the program and examples
of its usage and accuracy:

* L. Laaksonen, P. Pyykkö, and D. Sundholm, Fully numerical Hartree-Fock
  methods for molecules, Comp. Phys. Reports 4 (1986) 313-344. 
  http://doi.org/10.1016/0167-7977(86)90021-3

* L. Laaksonen, D. Sundholm, and P. Pyykkö, in "Scientific Computing in Finland",
  Eds. K. Kankaala and R. Nieminen, Research Report R1/89, Centre for Scientific Computing,
  (1989) p. 183.

* P. Pyykkö, in Numerical Determination of the Electronic Structure of Atoms, Diatomic
  and Polyatomic Molecules (NATO ASI Series C271) Eds. M. Defranceschi and J. Delhalle,
  (1989) p. 161.  http://doi.org/10.1007/978-94-009-2329-4

* J. Kobus, Finite-difference versus finite-element methods, Chem. Phys. Lett. 202
  (1993) 7-12.  http://doi.org/10.1016/0009-2614(93)85342-L

* J. Kobus, Vectorizable algorithm for the (multicolour) successive overrelaxation method,
   Comput. Phys. Commun. 78 (1994) 247-255. http://doi.org/10.1016/0010-4655(94)90003-5

* J. Kobus, L. Laaksonen, D. Sundholm, A numerical Hartree-Fock program for diatomic
   molecules, Comp. Phys. Commun. 98 (1996) 346-358. http://doi.org/10.1016/0010-4655(96)00098-7

* J.Kobus, D.Moncrieff, and S.Wilson, A comparison of the electric moments obtained from
  finite basis set and finite difference {Hartree-Fock} calculations for diatomic
  molecules, Phys. Rev. A 62 (2000), 062503/1--9

* J.Kobus, D.Moncrieff, and S.Wilson, Comparison of the polarizabilities and
  hyperpolarizabilities obtained from finite basis set and finite field, finite
  difference {Hartree-Fock} calculations for diatomic molecules, J. Phys. B:
  At. Mol. Opt. Phys. 34 (2001) 5127-5143.

* J. Kobus, Numerical Hartree-Fock methods for diatomic molecules, Handbook of Molecular
  Physics and Quantum Chemistry (Chichester), ed. S. Wilson (Wiley, 2002)

* J. Kobus, Hartree-Fock limit values of multipole moments, polarizabilities and
  hyperpolarizabilities for atoms and diatomic molecules, Comp. Lett. 3 (2007) 71-113.
  http://doi.org/10.1163/157404007782913408

* J. Kobus, A finite difference Hartree-Fock program for atoms and diatomic molecules,
  Comp. Phys. Commun. 184 (2013) 799-811. http://dx.doi.org/10.1016/j.cpc.2012.09.033

* J.Kobus, Hartree-Fock limit values of multipole moments, polarizabilities and
  hyperpolarizabilities for atoms and diatomic molecules, Phys. Rev. A 91 (2015) 022501,
  URL: http://link.aps.org/doi/10.1103/PhysRevA.91.022501,
  DOI: https://doi.org/10.1103/PhysRevA.91.022501

Density functionals are evaluated with the libxc library, which has been described in

* S. Lehtola et al, SoftwareX 2018, 7, 1–5.

The programming language used is Fortran 90 and C (if the multi-threaded version
employing pthreads is requested). Therefore Fortran 90 and C compilers together with
CMake are needed to compile and build the program. This can be done by running
`./x2dhfctl -b`.  See `./x2dhfctl -h` and `INSTALL` for more details.

`./doc/users-guide.pdf` contains the description of the program's input data and its
usage.

`test-sets/` directory contains the dozens of specific examples and the testctl
script should be used to list and run them.

`./lda_orbitals` and `./hf_orbitals` directories contain LDA and HF orbitals for a
number of atomic systems. These orbitals can be used to start the SCF process. The data
were generated by means of the HelFEM program and a finite-difference HF program for
atoms (qrhf), respectively,

`./bin` directory contains the x2dhf executable(s) and several BASH and Perl scripts
to facilitate the usage of the program (run `source .x2dhfrc` to adjust the
PATH variable accordingly). Firts of all try `xhf -h`. See also `elpropctl -h`,
`lxcctl -h`, `pecctl -h` and `testctl -h`.

To give the program a try, execute the following command: `cd tests; testctl run h/set-01.`

Jacek Kobus

2023-08-25

