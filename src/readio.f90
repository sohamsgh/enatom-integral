MODULE readio

    USE declarations

    CONTAINS

        subroutine READINPUT(c1, c2)
                ! read input file
                implicit none

                real(DP), intent(out)           :: c1, c2
                character(LEN=256)              :: infile

                infile = trim("avg.in")
                open(61, file=infile, status='old')
                read(61, '(F11.5)') c1
                read(61, '(F11.5)') c2
                close(61)

        end subroutine

        function READENTATOM(c1,c2) result (f)
                ! read the G-space function from enatom and 
                ! fourier transform it into real space

                ! for starters, I will use some analytic function
                ! This function is desgined to mimic the potential energy 
                ! of a single H sitting in the middle of the cell
                implicit none

                real(DP), intent(in)                    :: c1, c2

                integer                                 :: ix, iy, iz
                real(DP)                                :: r
                complex(DP),dimension(nr1,nr2,nr3)      :: f
                
                do iz = 1, nr3
                        do iy = 1, nr2
                                do ix = 1, nr1
                                        r = SQRT((ix*(a1-1)/nr1)**2 + (iy*(a2-1)/nr2)**2 + (iz*(a3-1)/nr3)**2)
                                        f(ix,iy,iz) = ((r-c1)/c2)*(1/(c2*SQRT(2*PI)))*EXP(-0.5*(((r - c1)/c2)**2))
                                end do
                        end do
                end do

        end function


        subroutine READEVAL(del, index_array, e_array)
              
        ! Reads the xml files in the K***** subdirectories of the prefix.save directory
        ! to find out bands within "del" of Fermi energy and their energy.
        ! Since the bands are arranged in ascending order of energy, 
        ! we can return the lowest bands index and the highest band index that satisfy the criterion.
        ! The information is arranged in a array which has nks columns and 2 rows. 
        ! The columns stand for each k-point. 
        ! If for a certain k-point, no band satisfies the criterion, then both rows are give -1. 
        ! If for a certain band, one or more bands satisfy the criterion, then the first row contains the 
        ! index of the lowest such band and the 2nd row contains the index of the highest such band
        ! Also calculates the energy of each band within the energy window. In this case no such trick is possible.
        ! This array is (nks, nbnd) in size. 
              implicit none

              real(DP), intent(in)                      :: del                          ! energy window around Fermi level
              integer, dimension(2,nks), intent(out)    :: index_array  ! stores the indices of k and bands that satisfy
              real(DP),dimension(nbnd,nks),intent(out)  :: e_array      ! eigenvalues of k and bands
              real(DP)                                  :: eval                         ! eigenvalue
              character(len=256), parameter             :: catchphrase1 = "<EIGENVALUES"
              character(len=256), parameter             :: catchphrase2 = "</EIGENVALUES>"
              character(len=256)                        :: fname                        ! filename containing the eigenvalue xml data
              character(len=256)                        :: temp                         ! temp reads string headers in xml file
              integer                                   :: inks,inbnd,buffer            ! iterators for k-point and bands
              character(len=8)                          :: fmt                          ! format descriptor
              character(len=5)                          :: x1                           ! for creating the write file name for each kpoint      
              logical                                   :: state, flag
              fmt = '(I5.5)'
              
              !write(6, '(A, 4X, F7.4)') "Fermi Energy in Rydberg =", EFermi
              !write(6, '(A)') "Reading the eigenvalue xml files..."
              !write(6, '(A)') "----------------------------------------"
              index_array = (-1)                      ! Initial array to be a value that will be rejected on reading
              e_array = (1000000.00)                       ! some absurdly high number
              do inks = 1, nks
                write(x1, fmt) inks             ! converting integer to string using an 'internal file'
                fname = trim(tmpdir)//trim(prefix)//".save/"//"K"//trim(x1)//"/"//"eigenval.xml"
                state = .TRUE.
                buffer = 0
        
                open(601, file=fname, status='old')
                !write(6, '(A, A)') "Reading file", trim(fname)  
        
                do while (state .eqv. .TRUE.)
                        
                        read(601, 990) temp
                        990 format (A60)
                        temp = trim(adjustl(temp))
                        buffer = buffer + 1
                        if (buffer .ge. 100000) then
                                !write(6, '(A, A)') "infinite loop encountered while reading from file  = ",fname
                                !write(6, '(A, A)') "xml catchphrase to end reading not found", catchphrase2
                                STOP 0
                        end if
        
                        if (temp(1:12) .eq. "<EIGENVALUES") then
                                ! Locate Eigenvalue block
                                flag = .TRUE.
                                do inbnd = 1, nbnd
                                        read(601, *) eval
                                        e_array(inbnd, inks) = eval
                                        if (ABS(eval - EFermi) .LE. del) then
                                                if (flag .eqv. .TRUE.) then    ! This should occur only once
                                                        !write(6, '(A, 4X, F10.6, 4X, A, 4X, I5)') &
                                                                !&"eval in =", eval, "lower band =", inbnd
                                                        ! lowest band that falls within energy window
                                                        index_array(1,inks) = inbnd
                                                        index_array(2,inks) = inbnd   ! set highest band the same as lowest
                                                        flag = .FALSE.   ! The only place flag is set to .FALSE. Should be so.
                                                        ! The following if condition is executed only once, for flag = TRUE   
                                                        ! the equality condition should hold if the flag and index_array assignments
                                                        ! are done correctly.
                                                        if (e_array(index_array(1,inks), inks) .ne. eval) then
                                                                !write(6, '(A)') "Something wrong with energies"
                                                                STOP 0
                                                        end if
                                                else                            ! This should should occur m-1 times, where
                                                                                ! there are m bands in the window.
                                                        index_array(2, inks) = inbnd
                                                end if
                                        end if
                                end do
                                if (index_array(2,inks) .gt.0) then
                                        !write(6, '(A, 4X, F10.6, 4X, A, 4X, I5)') "eval_out &
                                                !&=", e_array(index_array(2, inks), inks), &
                                                !&"higher band = ", index_array(2,inks)
                                end if
                                state = .FALSE.
                        !write(6, '(A)') "----------------------------------------"
                        end if

                end do
        
                close(601)
              end do
        end subroutine


END MODULE readio
