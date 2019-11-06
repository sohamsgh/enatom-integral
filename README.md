# enatom-integral
Fortran code for calculating double integrals with a (k,k') kernel on the Fermi Surface.

The computational scheme has two braod parts.
  1. Creating the kernel, which is a function of k and k' vectors. This is done through the enatom procedure. 
  2. Solving the double integral on the Fermi Surface with this kernel, using terahedral interpolation.
  
  declarations.f90:  Declarations of the constants of the problem.
  
  compute.f90:       Computational Subroutines.
  
  readio.f90:        input/output procedures.
  
  dos.f90:           calculate DOS using tetrahedron integration.
  
