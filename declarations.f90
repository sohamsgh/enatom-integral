      module declarations

      ! Store constants of the system being calculated

      implicit none

      save

      ! DP
      integer DP
      parameter (DP=kind(1.d0))

      ! DIRECT_TO_FACTOR = a binary byte factor which is
      ! Different between ifort and gfortran
      integer DIRECT_IO_FACTOR
      parameter (DIRECT_IO_FACTOR=8)

      ! PI
      real(DP) PI
      parameter (PI=4.D0*DATAN(1.D0))

      ! k-mesh, shift in BZ origin 
      integer n1, n2, n3, orig1, orig2, orig3
      parameter (n1=20,n2=20,n3=20, orig1=0, orig2=0, orig3=0)

      !number of k points, bands, dimensions.
      integer nks, nbnd, ic
      parameter (nks=n1*n2*n3,nbnd=4,ic=3)

      !real space lattice constants (Bohr)
      real(DP) a1, a2, a3
      parameter (a1=8.6666667,a2=8.6666667,a3=8.6666667)

      !Fermi energy (Hartree)
      real(DP) EFermi
      parameter (EFermi=0.0)
      !parameter (EFermi=-2.309547862172711E-001)

      !real space wavefunction file
      character(len=256) wfcrfile
      parameter(wfcrfile="h_jellium.wfc_r")

      !temp directory, prefix
      character(len=256) tmpdir, prefix
      parameter (tmpdir="./", prefix="h_jellium")

      !real space grid (given as FFT grid in scf output).
      integer nr1, nr2, nr3
      parameter (nr1=45,nr2=45,nr3=45)
    
      !width of the energy window around Fermi level (eV)
      real(DP) del
      parameter (del=5.0)

      
      end module declarations

