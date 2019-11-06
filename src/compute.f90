MODULE compute  

      use declarations
      use FFTW3
      USE OMP_LIB

      implicit none

      CONTAINS

        subroutine index3dto1d(id, n2, n3, id1, id2, id3)
        ! Given id and (n1, n2, n3) return id1, id2, id3. 
        ! n1, n2, n3 are range of id1, id2, id3. id1 ranges from 1 to n1 etc.
        ! id3 is taken to move the fastest, id1 the slowest.


                integer, intent(in)             :: id
                integer, intent(in)             :: n2, n3
                integer, intent(out)            :: id1, id2, id3
            

                ! from given index, calculate x-index, y-index, z-index
                id1 = INT((id-1)/(n2*n3))
                id2 = INT((id - 1 - id1*n2*n2)/n3)
                id3 = MOD(id-1,n3)

                id1 = id1 + 1
                id2 = id2 + 1
                id3 = id3 + 1

        end subroutine

        pure function index1dto3d(id1, id2, id3, n2, n3) result(id)
        ! Given id and (n1, n2, n3) return id1, id2, id3. 
        ! n1, n2, n3 are range of id1, id2, id3. id1 ranges from 1 to n1 etc.
        ! id3 is taken to move the fastest, id1 the slowest.

                integer, intent(in)             :: id1, id2, id3
                integer, intent(in)             :: n2, n3
                integer                         :: id
            

                ! from given x-index, y-index, z-index, calculate index

                id = (id1-1)*n2*n3 + (id2-1)*n3 + id3

        end function


        function indextok(id, n1, n2, n3) result(k)
        ! given an index of a 3D vector k, gives the vector k as a (3) array 
        ! example 1. if k = [0/7, 1/7, 2/7, 3/7, -3/7, -2/7, -1/7]

                integer, intent(in)             :: id
                integer, intent(in)             :: n1, n2, n3
                
                integer, dimension(3)           :: id1d
                real(DP), dimension(3)          :: k
                integer, dimension(3)           :: n
                integer                         :: i


                call index3dto1d(id, n2, n3, id1d(1), id1d(2), id1d(3))
                n(1) = n1; n(2) = n2; n(3) = n3

                do i = 1, 3
                        if (id1d(i) .le. INT((n(i)+1)/2)) then
                                k(i) = REAL(id1d(i)-1)/n(i)
                        else if (id1d(i) .gt. INT((n(i)+1)/2)) then
                                k(i) = REAL(id1d(i)-1)/n(i) -1
                        else
                                print*, "In function indextok"
                                print*,"There is something wrong with the value of k-index:", id1d(i)
                        STOP 0
                        end if
                end do

        end function

        function ktoindex(k, n1, n2, n3) result(id)
        ! given a 3D vector k, gives its index between 1 and n
        ! example 1. if k_1 = [0/8, 1/8, 2/8, 3/8, -4/8, -3/8, -2/8, -1/8]
        !              id_1 = [1, 2, 3, 4, 5, 6, 7, 8]
        ! and similarly for k_2 --> id_2 and k_3 --> id_3. 
        ! then from id_1, id_2, id_3, id is created
        ! Note the use of NINT - nearest integer, rather than INT.

                real(DP), dimension(3), intent(in)      :: k
                integer, intent(in)                     :: n1, n2, n3
                integer                                 :: id

                integer, dimension(3)                   :: id1d
                integer, dimension(3)                   :: n
                integer                                 :: i
                n(1) = n1; n(2) = n2; n(3) = n3
                !if (n*k .NE. INT(n*k)) then 
                !        print*, "In function indextok"
                !        print*, "n*k = ", n*k, "is not an integer, that sounds wrong"
                !        STOP 0
                !end if
                do i = 1, 3
                        if (int(n(i)*k(i)) .ge. 0) then
                                id1d(i) = NINT(n(i)*k(i)) + 1
                        else
                                id1d(i) = NINT(n(i)*(k(i)+1)) + 1
                                if (id1d(i) .GT. n(i) ) then
                                        write(6,'(A, 2X, I4, F8.3, I8)') "problem:", i, k(i), id1d(i)
                                end if
                        end if
                end do
                !write(6, *) id1d(1), id1d(2), id1d(3)
                id = index1dto3d(id1d(1), id1d(2), id1d(3), n2, n3)

        end function
        
        subroutine scan_BZ(k, dk, c_m)
        ! Given a k-point inside the BZ and the grid spacing, 
        ! create the corners of the submesh cell.
        ! k is designated "vertex 3" in accordance with Blochl 1994.
        ! And then move in the +ve x, y, z direction one step
        ! Question: Why do we need both k and (ik1, ik2, ik3)?
                implicit none
                
                real(DP), intent(in)                    :: k(:), dk(:)
                real(DP), intent(out)                   :: c_m(:,:)

                integer                                 :: i, ik1, ik2, ik3
                real(DP)                                :: dk1, dk2, dk3
        
                ik1 = NINT(k(1)*n1)
                ik2 = NINT(k(2)*n2)
                ik3 = NINT(k(3)*n3)
                ! Consider the special cases of the boundary
                if ( ik1 .EQ. (n1-1)/2 ) then 
                        dk1 = -REAL(n1-1)/n1
                else
                        dk1 = dk(1)
                end if
                if ( ik2 .EQ. (n2-1)/2 ) then 
                        dk2 = -REAL(n2-1)/n2
                else
                        dk2 = dk(2)
                end if
                if ( ik3 .EQ. (n3-1)/2 ) then 
                        dk3 = -REAL(n3-1)/n3
                else
                        dk3 = dk(3)
                end if
                ! Create the corners of the parallelopiped submesh
                c_m(:,3) = k
                c_m(:,1) = (/ k(1) + dk1, k(2), k(3) /)
                c_m(:,4) = (/ k(1), k(2) + dk2, k(3) /)
                c_m(:,7) = (/ k(1), k(2), k(3) + dk3 /)
                c_m(:,2) = (/ k(1) + dk1, k(2) + dk2, k(3) /)
                c_m(:,8) = (/ k(1), k(2) + dk2, k(3) + dk3 /)
                c_m(:,5) = (/ k(1) + dk1, k(2), k(3) + dk3 /)
                c_m(:,6) = (/ k(1) + dk1, k(2) + dk2, k(3) + dk3 /)

        end subroutine 
       

        subroutine scan_BZ2(ik1, ik2, ik3, k, dk, c_m)
        ! Given a k-point inside the BZ and the grid spacing, 
        ! create the corners of the submesh cell.
        ! k is designated "vertex 3" in accordance with Blochl 1994.
        ! And then move in the +ve x, y, z direction one step
        ! Question: Why do we need both k and (ik1, ik2, ik3)?
                implicit none
                
                integer, intent(in)                     :: ik1, ik2, ik3
                real(DP), intent(in)                    :: k(:), dk(:)
                real(DP), intent(out)                   :: c_m(:,:)

                real(DP)                                :: dk1, dk2, dk3
        
                ! Consider the special cases of the boundary
                if ( ik1 .EQ. (n1-1)/2 ) then 
                        dk1 = -REAL(n1-1)/n1
                else
                        dk1 = dk(1)
                end if
                if ( ik2 .EQ. (n2-1)/2 ) then 
                        dk2 = -REAL(n2-1)/n2
                else
                        dk2 = dk(2)
                end if
                if ( ik3 .EQ. (n3-1)/2 ) then 
                        dk3 = -REAL(n3-1)/n3
                else
                        dk3 = dk(3)
                end if
                ! Create the corners of the parallelopiped submesh
                c_m(:,3) = k
                c_m(:,1) = (/ k(1) + dk1, k(2), k(3) /)
                c_m(:,4) = (/ k(1), k(2) + dk2, k(3) /)
                c_m(:,7) = (/ k(1), k(2), k(3) + dk3 /)
                c_m(:,2) = (/ k(1) + dk1, k(2) + dk2, k(3) /)
                c_m(:,8) = (/ k(1), k(2) + dk2, k(3) + dk3 /)
                c_m(:,5) = (/ k(1) + dk1, k(2), k(3) + dk3 /)
                c_m(:,6) = (/ k(1) + dk1, k(2) + dk2, k(3) + dk3 /)

        end subroutine 

        function tetra(array_k, array_energy, array_quantity, E) result(I)
        ! given the indices of the four corners of a tetrahedron,
        ! and the values of a quantity on those corresponding corners, 
        ! calculate the integral inside the tetrahedron.
        !
        ! The order of the input indices are immaterial, but
        ! all the inputs should have the same order
        ! Limited to real functions.
        ! array_k is a 2D (3X4) array. Rows for 3 components k1, k2, k3. Columns are for 4 vertices
        ! array_energy and array_A are 1-D arrays of size (4).
        ! E is the constant energy (usually Fermi energy)

                implicit none

                real(DP), intent(IN)                :: array_k(:,:)
                real(DP), intent(IN)                :: array_energy(:)
                real(DP), intent(IN)                :: array_quantity(:)
                real(DP), intent(IN)                :: E
                real(DP)                            :: I
               
                real(DP), dimension(3)                                  :: n
                integer, allocatable, dimension(:)                      :: indices
                real(DP), allocatable, dimension(:,:)                   :: sorted_k, dummy_k         
                real(DP), allocatable, dimension(:)                     :: sorted_energy
                real(DP), allocatable, dimension(:)                     :: sorted_quantity
                real(DP), allocatable, dimension(:)                     :: a,b
                real(DP), allocatable, dimension(:,:)                   :: sorted_r
                real(DP)                                                :: vol, rspace_vol   

                real(DP)                                                :: f0binv, f1binv, f3binv
                real(DP), allocatable, dimension(:)                     :: s0, s1, s3, sbar
                real(DP)                                                :: i0
                real(DP), allocatable, dimension(:)                     :: i1
                real(DP), allocatable, dimension(:)                     :: k10, k20, k30, k21, k31, k32
                real(DP)                                                :: e10, e20, e30, e21, e31, e32
                real(DP)                                                :: ee0, ee1, ee2, ee3
                integer                                                 :: size1, size2
                integer                                                 :: iv, iloc
                real(DP)                                                :: tolerance

                size1 = size(array_k, 1); size2 = size(array_k, 2)  ! should be size1 = 3, size2 = 4
                if ((size1 .NE. 3) .AND. (size2 .NE. 4)) then
                        write(6,*) "Error in dimensions of input array array_k", array_k
                        STOP 1
                end if
                allocate(sorted_k(size1, size2), dummy_k(size1, size2))
                allocate(sorted_energy(size2))
                allocate(sorted_quantity(size2))
                allocate(indices(size2))
                allocate(a(size1), b(size1))
                allocate(sorted_r(size1, size2))
                allocate(i1(size1))
                allocate(s0(size1), s1(size1), s3(size1), sbar(size1))
                allocate(k10(size1),k20(size1),k30(size1),k21(size1),k31(size1),k32(size1)) 

                ! Initialize tetrahedron integral to zero
                I = 0.0
                ! initialize/set other quantities
                !
                n(1) = n1; n(2) = n2; n(3) = n3
                f0binv = 0.0; f1binv = 0.0; f3binv = 0.0
                s0 = 0.0; s1 = 0.0; s3 = 0.0
                i0 = 0 
                i1 = 0
                tolerance = 1.D0/10000000.0
                do iv = 1,size2
                        indices(iv) = iv
                end do
                ! copy input k-array to in-function array
                ! sorted_k = array_k
                ! copy input E-array to in-function array
                 sorted_energy = array_energy
                ! sort the energies in ascending order and get the order of indices so the k and A values can be ordered   
                call quicksort(sorted_energy, indices)

                do iv = 1, size2
                        dummy_k(:,iv) = array_k(:,indices(iv))
                        sorted_quantity(iv) = array_quantity(indices(iv))
                end do
                ! Make k0 the origin and other k-vectors w.r.t. it. 
                do iv = 2,size2
                        sorted_k(:,iv) = dummy_k(:,iv) - dummy_k(:,1)
                end do
                sorted_k(:,1) = 0.0

                ! Check for edge of the BZ
                do iv = 1, 4
                        do iloc = 1,3
                                if (sorted_k(iloc, iv) .GT. REAL(n(iloc)-2)/n(iloc)) then
                                        sorted_k(iloc,iv) = sorted_k(iloc,iv) - &
                                                & REAL(n(iloc) -2)/n(iloc)
                                else if (sorted_k(iloc, iv) .LT. -REAL(n(iloc)-2)/n(iloc)) then
                                        sorted_k(iloc,iv) = sorted_k(iloc,iv) + &
                                                & REAL(n(iloc) -2)/n(iloc)
                                else
                                end if
                        end do
                end do

                ! Calculate volume v = 6 times the volume of the tetrahedron.
                vol =  ABS(6.0*DOT_PRODUCT(sorted_k(:,2), cross(sorted_k(:,3), sorted_k(:,4))))
                ! The following condition should never be satisfied
                if (vol .LT. tolerance) then
                        write(6,*) "the tetrahedron has zero volume"
                        write(6,'(A, 3F8.4,2X,3F8.4,2X,3F8.4,2X,3F8.4)') "vertices:", sorted_k
                        STOP 1
                end if
                ! calculate the real space contragradients to k-vectors
                sorted_r(:,1) = 0.0   ! probably not even necessary to calculate this
                sorted_r(:,2) = cross(sorted_k(:,3), sorted_k(:,4))/vol
                sorted_r(:,3) = cross(sorted_k(:,4), sorted_k(:,2))/vol
                sorted_r(:,4) = cross(sorted_k(:,2), sorted_k(:,3))/vol
                rspace_vol =  ABS(6.0*DOT_PRODUCT(sorted_r(:,2), cross(sorted_r(:,3), sorted_r(:,4))))
                if (rspace_vol .LT. tolerance) then
                        write(6,*) "the tetrahedron has zero real space volume"
                        write(6,*) "vertices:", sorted_r
                        STOP 1
                end if
                  
                ! calculate a and b
                a = 0.D0
                b = 0.D0
                ! calculate quantities in the denominator

                do iv = 2, size2
                        b(:) = b(:) + (sorted_energy(iv) - sorted_energy(1))*sorted_r(:,iv)
                        a(:) = a(:) + (sorted_quantity(iv) - sorted_quantity(1))*sorted_r(:,iv)
                end do


                ! try the different possible locations of energies
                ! Note the If statements designed to avoid zeros in the denominators. 
                ! Let's calculate all the quantities first 
                ! and then go onto the if conditions. Prioritise readability over slight speed 
                ! improvement. 
                ! Also sticking to Lehman and Taut's monenclature, that means one of two superflous assignments
                

                e10 = sorted_energy(2) - sorted_energy(1)
                e20 = sorted_energy(3) - sorted_energy(1)
                e30 = sorted_energy(4) - sorted_energy(1)
                e21 = sorted_energy(3) - sorted_energy(2)
                e31 = sorted_energy(4) - sorted_energy(2)
                e32 = sorted_energy(4) - sorted_energy(3)
                ee0 = E - sorted_energy(1)
                ee1 = E - sorted_energy(2)
                ee2 = E - sorted_energy(3)
                ee3 = E - sorted_energy(4)
                ! calculate other quantity related stuff
                k10 = sorted_k(:,2) - sorted_k(:,1)
                k20 = sorted_k(:,3) - sorted_k(:,1)
                k30 = sorted_k(:,4) - sorted_k(:,1)
                k21 = sorted_k(:,3) - sorted_k(:,2)
                k31 = sorted_k(:,4) - sorted_k(:,2)
                k32 = sorted_k(:,4) - sorted_k(:,3)

                if  ((sorted_energy(1) < E) .AND. (E <= sorted_energy(2))) then
                        ! if e10 becomes very small, ee0 is also small since it falls in between.
                        ! (in the case E ~ e1 ~ e2). In that case set the surface values to zero.
                        if (ABS(e10) .LT. tolerance) then
                                ! calculate DOS related stuff
                                f0binv = 0.D0
                                ! calculate other quantity related stuff
                                s0 = sorted_k(:,1)
                        else
                                ! calculate DOS related stuff
                                f0binv = (vol/2.D0)*((ee0**2)/(e10*e20*e30))
                                ! calculate other quantity related stuff
                                s0 = sorted_k(:,1) + (ee0/3.D0)*(k10/e10 + k20/e20 + k30/e30)
                        end if
                        i0 = f0binv      ! reassign variable to make them the same as Lehmann and Taut
                        i1 = f0binv*s0
                        ! the following condition is tricky. Check notes. 
                else if  ((sorted_energy(2) < E) .AND. (E <= sorted_energy(3))) then
                        if (ABS(e10) .LT. tolerance) then
                                f0binv = 0.D0
                                f1binv = 0.D0
                                s0 = sorted_k(:,1)
                                s3 = sorted_k(:,4)
                                sbar = 0.D0
                        else if (ABS(e21) .LT. tolerance) then
                                f0binv = (vol/2.D0)*((ee0**2)/(e10*e20*e30))
                                f1binv = 0.D0
                                sbar = 0.D0
                        else
                                f0binv = (vol/2.D0)*((ee0**2)/(e10*e20*e30))
                                f1binv = (vol/2.D0)*((ee1**2)/(e10*e21*e31))
                                ! k01/e01 is equal to k10/e10
                                s0 = sorted_k(:,1) + (e10/3.D0)*(k10/e10 + k20/e20 + k30/e30) 
                                s3 = sorted_k(:,4) + (-e32/3.D0)*(k30/e30 + k31/e31 + k32/e32) 
                                sbar = (-ee2/e21)*s0 + (ee1/e21)*s3
                        end if
                        i0 = f0binv -f1binv
                        i1 = i0*sbar
                                
                        !i1 = 0.D0
                else if  ((sorted_energy(3) <= E) .AND. (E < sorted_energy(4))) then
                        if (ABS(e32) .LT. tolerance) then
                                f3binv = 0.D0
                                s3 = 0.D0
                        else
                                f3binv = (vol/2.0)*((ee3**2)/(e30*e31*e32))
                                ! k31/e31 equal k13/e13 etc. 
                                s3 = sorted_k(:,4) + (ee3/3.0)*(k30/e30 + k31/e31 + k32/e32)
                        end if
                        i0 = f3binv                                ! another such harmless step
                        i1 = f3binv*s3
                else
                        !write(6,'(A)') "Condition 4"
                        i0 = 0
                        i1 = 0
                end if

                I = sorted_quantity(1)*i0 + (a(1)*i1(1) + a(2)*i1(2) + a(3)*i1(3))
                !I = sorted_quantity(1)*i0 
                
                deallocate(indices)
                deallocate(sorted_k, dummy_k)
                deallocate(sorted_energy)
                deallocate(sorted_quantity)
                deallocate(a, b)
                deallocate(sorted_r)
                deallocate(i1)
                deallocate(s0, s1, s3)
                deallocate(k10,k20,k30,k21,k31,k32) 



        end function tetra

        function cross(a, b) result(c)

                implicit  none
        
                real(DP), intent(in)            :: a(:), b(:)
                real(DP)                        :: c(3)

                if ((size(a) .ne. 3) .OR. (size(b)) .ne. 3) then
                        print*,"cross product only for three dimensional vectors"
                        print*,"rank of vector1:", size(a), "rank of vector 2:", size(b)
                        stop 1
                end if

                c(1) = a(2)*b(3) - a(3)*b(2)
                c(2) = a(3)*b(1) - a(1)*b(3)
                c(3) = a(1)*b(2) - a(2)*b(1)
        
        end function cross

        pure function trapezoid(x,y) result(r)
        !! Calculates the integral of an array y with respect to x using the trapezoid approximation. 
        !! The mesh x does not have to be uniform.
        real(DP), intent(in)                    :: x(:)         !! variable x
        complex(DP), intent(in)                 :: y(size(x))   !! variable y, complex
        complex(DP)                             :: r            !! Integral y . dx

                ! Integrate using the trapezoid approximation
                associate (n => size(x))
                        r = sum((y(1+0:n-1) + y(1+1:n-0))*(x(1+1:n-0) - x(1+0:n-1)))/2
                end associate
        end function

        function trapz_parallel(x,y) result(r)
        !! Calculates the integral of an array y with respect to x using the trapezoid approximation. 
        !! The mesh x does not have to be uniform.
        !! parallelizes via OpenMP.
        real(DP), intent(in)                    :: x(:)         !! variable x
        complex(DP), intent(in)                 :: y(size(x))   !! variable y, complex
        complex(DP)                             :: r            !! Integral y . dx
        integer                                 :: i, N         !! index, size of array

        r = (0.D0, 0.D0)
        N = size(x)
        !$omp parallel do shared(x,y) reduction(+:r)
            do i = 1, N-1
                r = r + ((y(i) + y(i+1))*(x(i+1) - x(i)))/2
            end do
        !$omp end parallel do

        end function


        recursive subroutine quicksort(a, indices)
        ! What it says. Its a increasing order quicksort code.
        ! Also produces an array of indices ordered according to the sorted elements
        ! so if a = [3,1,4,2], then it has by default indices = [1,2,3,4]. After sorting, 
        ! indices = [2, 4, 1, 3]. But indices could be any other integer array as well.

                implicit none

                real(DP)                :: a(:)
                integer, optional       :: indices(:)
                real(DP)                :: x, t
                integer                 :: m
                integer                 :: first = 1, last
                integer                 :: i, j

                last = size(a, 1)
                if (size(indices,1) .ne. last) then
                        print*, " Size mismatch between array and indices"
                        stop 1
                end if

                x = a((first+last)/2)
                i = first
                j = last

                do 
                        do while (a(i) < x)
                                i = i + 1
                        end do
                        do while (x < a(j))
                                j = j - 1
                        end do
                        if (i >=j) exit
                        t = a(i); a(i) = a(j); a(j) = t
                        if (present(indices)) then
                                m = indices(i); indices(i) = indices(j); indices(j) = m
                        end if
                        i = i + 1
                        j = j -1
                end do

                if (first < i-1) call quicksort(a(first:i-1), indices(first:i-1))
                if (j+1 < last)  call quicksort(a(j+1:last), indices(j+1:last))
        
        end subroutine quicksort
        
        subroutine chk_tetra(k1,k2,k3,k4,e_array,min_band,max_band,array_tetra_index,array_tetra,E,loci)
        ! takes in as input 4 k-points, the (being created) array 
        ! which holds the references of all the tetrahedrons that are active, as in have the Fermi 
        ! surface (or any other constant energy surface going through it, and an array 
        ! that holds the corners and the energies themselves. The two  arrays 
        ! has as many columns as the number of active tetrahedrons. In each column, 
        ! the index array has the id of the four corners, then one band index each of the 
        ! four k-point corners.
        ! The array_tetra holds the corresponding real numbers. The k-vectors themselves, 
        ! and then the four energy values at the corners corresponding to those bands. 
        ! min_band and max_band are the lowest and highest bands that cross the Fermi surface.

        ! loci is the column number of the arrays where the data is to be enetered. 
        ! loci is an 'inout' quantity. Every time there is a tetrahedron that is cut by the constant 
        ! energy surface, loci is increased by one. We are considering the same four k-point corners 
        ! but different bands, hence different energies at the corners, as different tetrahedrons. 
        ! So at the end of all the calls to this subroutines, loci should be the number of total 
        ! active tetrahedrons. 

                implicit none
               
                real(DP), dimension(:), intent(in)      :: k1, k2, k3, k4
                real(DP), intent(in)                    :: e_array(:,:)
                integer, intent(in)                     :: min_band, max_band
                integer, intent(inout)                  :: array_tetra_index(:,:)
                real(DP), intent(inout)                 :: array_tetra(:,:)
                integer, intent(inout)                  :: loci  ! counter for num of columns (=num of active tetras)

                integer                                 :: id1, id2, id3, id4
                integer                                 :: ibnd1, ibnd2, ibnd3, ibnd4
                real(DP)                                :: iemin, iemax
                real(DP)                                :: E, e21, e31, e41
                real(DP), dimension(4)                  :: e_vertex
                real(DP)                                :: tolerance

                tolerance = 1.D0/1000000000.0

                ! Check 1: A quick check to count out a class of tetrahedrons
        
                ! If the lowest band of each of the 4 k-points is above the Fermi 
                ! level or if the highest band of each of the k-points is below 
                ! the Fermi level, that tetrahedron (those 4 k-points, all bands) does not 
                ! contribute.
                id1 = ktoindex(k1, n1, n2, n3)
                id2 = ktoindex(k2, n1, n2, n3)
                id3 = ktoindex(k3, n1, n2, n3)
                id4 = ktoindex(k4, n1, n2, n3)
                iemin = min(e_array(min_band,id1),e_array(min_band,id2),e_array(min_band,id3),e_array(min_band,id4))  
                iemax = max(e_array(max_band,id1),e_array(max_band,id2),e_array(max_band,id3),e_array(max_band,id4))  
                if ( (iemin .le. E) .AND. (iemax .ge. E) ) then   ! this means the tetrahedron passes Check 1
                        ! now we have to check individual bands
                        ! Get the k-vectors of the corners from their indices
                        ! Check 2. Verify energy of each band individually
                        do ibnd1 = min_band, max_band
                            do ibnd2 = min_band, max_band
                                do ibnd3 = min_band, max_band
                                    do ibnd4 = min_band, max_band
                                        e_vertex(1) = e_array(ibnd1, id1)
                                        e_vertex(2) = e_array(ibnd2, id2)
                                        e_vertex(3) = e_array(ibnd3, id3)
                                        e_vertex(4) = e_array(ibnd4, id4)
                                        e21 = e_vertex(2) - e_vertex(1)
                                        e31 = e_vertex(3) - e_vertex(1)
                                        e41 = e_vertex(4) - e_vertex(1)
                                        ! Check 2. On whether the tetrahedron should be included or not
                                        if ( (minval(e_vertex) .le. E) .AND. (maxval(e_vertex) .ge. E) ) then 
                                            ! There is no apriori way to know how many times this if condition 
                                            ! but every time it does, we have one new active tetra
                                            loci = loci + 1
                                            array_tetra_index(1:4, loci) = (/ id1, id2, id3, id4 /)
                                            array_tetra_index(5:8, loci) = (/ ibnd1, ibnd2, ibnd3, ibnd4 /)
                                            array_tetra(1:3, loci) = k1
                                            array_tetra(4:6, loci) = k2
                                            array_tetra(7:9, loci) = k3
                                            array_tetra(10:12, loci) = k4
                                            array_tetra(13:16, loci) = e_vertex
                                            !write(6,'(3I8, F8.5)' INT(n1*array_tetra(1:3, loci), & 
                                            !        & array_tetra(13, loci)
                                            !write(6, *) "Passed Check 1 and 2 for corners:", id1, id2, id3, id4
                                            !write(6,'(A, F6.4, A, F6.4, A, F6.4 )') &
                                            !        & "min(e_vertex)=  ", minval(e_vertex), "  Energy= ", E, &
                                            !        & "  max(e_vertex) = ", maxval(e_vertex)
                                            !write(6, *) "loci =", loci, "  array_tetra_index =", &
                                            !        & array_tetra_index(:,loci)
                                        else 
                                            !write(6, *) "I have passed check 1 but not &
                                            !        & check 2 for indices", id1, id2, id3, id4 

                                        end if
                                    end do
                                end do
                            end do
                        end do
                else
                        !write(6,*) "                         "
                        !write(6,*) "E is outside the tetrahedra for all bands"
                end if
        end subroutine chk_tetra


        function expect(dist_evc_r_ik, dist_evc_r_ikp, func) result (weight)
                ! return the expectation squared value of the function func(r) 
                ! using  integral dxdydz wfc_(k'n') (x,y,z) func(x,y,z) wfc_(k,n)(x,y,z)
                ! function func itself comes from enatom code.
                ! The wavefunctions are fiven in 1-D form, with 
                ! (assuming) x index increasing the fastest and z the slowest. 

                implicit none

                complex(DP), intent(in)         :: dist_evc_r_ik(:), dist_evc_r_ikp(:) ! 3D wfc in 1-D form
                complex(DP), intent(in)            :: func(:,:,:)  ! f(x,y,z)
                real(DP)                        :: vol, weight            ! squared expectation

                COMPLEX(DP), allocatable        :: wfc_block_ik(:), wfc_block_ikp(:) !psi(x) for fixed y and z
                COMPLEX(DP), ALLOCATABLE        :: M1(:), M2(:,:)
                REAL(DP), ALLOCATABLE           :: xarray(:), yarray(:), zarray(:)
                REAL(DP), ALLOCATABLE           :: sum_M1(:)
                integer                         :: xloc  ! location to get psi(:, iy, iz)                
                integer                         :: i, iy, iz  
                
                ALLOCATE(M2(nr2, nr3))
                ALLOCATE(wfc_block_ik(nr1), wfc_block_ikp(nr1), M1(nr3), sum_M1(nr3))
                ALLOCATE(xarray(nr1), yarray(nr2), zarray(nr3))

                xarray = (0.d0) 
                yarray = (0.d0) 
                zarray = (0.d0)
                vol = a1*a2*a3

                do i = 1, nr1
                        xarray(i) = a1*(i-1)/nr1 + 0.001
                end do
                do i = 1, nr2
                        yarray(i) = a2*(i-1)/nr2 + 0.001
                end do
                do i = 1, nr3
                        zarray(i) = a3*(i-1)/nr3 + 0.001
                end do

                wfc_block_ik=(0.d0,0.d0)
                wfc_block_ikp=(0.d0,0.d0)
                M2 = (0.0,0.0)
                M1 = (0.0,0.0)
                
                DO iz = 1,nr3
                        DO iy = 1,nr2
                                xloc = nr1*(iy-1) + nr1*nr2*(iz-1)
                                wfc_block_ik = dist_evc_r_ik((xloc+1): (xloc+nr1))/sqrt(vol)
                                wfc_block_ikp = dist_evc_r_ikp((xloc+1): (xloc+nr1))/sqrt(vol)
                                ! Call trapezoidal integration
                                M2(iy,iz)=trapezoid(xarray,CONJG(wfc_block_ikp)*func(:,iy, iz)*wfc_block_ik)
                        END DO
                        M1(iz) = trapezoid(yarray, M2(:,iz))
                END DO
                weight = (ABS(trapezoid(zarray, M1)))**2
                
                DEALLOCATE(wfc_block_ik, wfc_block_ikp, M1,M2, sum_M1)
                DEALLOCATE(xarray, yarray, zarray)

        end function 

        function fft3d(threedin, L, M, N, norm) result(threedout)
                ! implement fast fourier transform using the fftw3 module

                complex(C_DOUBLE_COMPLEX), intent(in)           :: threedin(:,:,:)
                integer, intent(in)                             :: L, M, N
                
                complex(C_DOUBLE_COMPLEX), dimension(L, M, N)   :: dummy_in, threedout
                
                type(C_PTR)                                     :: plan
                integer                                         :: i, j, k
                logical                                         :: norm
                integer                                         :: weight 
                
                dummy_in = threedin
                ! note the reversal of the indices L, M, N
                plan = fftw_plan_dft_3d(N, M, L, dummy_in,threedout, FFTW_FORWARD,FFTW_ESTIMATE)

                call fftw_execute_dft(plan, dummy_in, threedout)

                weight = L*M*N

                if ( norm .EQV. .TRUE.) then
                        ! normalize the output data
                        do k = 1, N
                                do j = 1, M
                                        do i = 1, L
                                                threedout(i,j,k) = threedout(i,j,k)/weight
                                        end do
                                end do
                        end do
                end if

                call fftw_destroy_plan(plan)

        end function



END MODULE compute
