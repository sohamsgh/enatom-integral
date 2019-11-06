PROGRAM INTEGRATION_CHECKER
        ! Check the trapezoidal integration function
        ! As well as its parallel version.
        use declarations
        use compute
        use OMP_LIB
        IMPLICIT NONE

        REAL(DP), ALLOCATABLE           :: array_x(:)
        COMPLEX(DP), ALLOCATABLE        :: array_y(:)
        INTEGER                         :: n, i
        REAL(DP)                        :: a, b, x, dx, r
        COMPLEX(DP)                     :: y


        n = 1000000000

        ALLOCATE(array_x(n))
        ALLOCATE(array_y(n))

        a = 0.D0; b = 1.D0
        dx = (b-a)/n
       !$omp parallel do shared(array_x,array_y) 
        do i = 0, n-1
                x = a + i*dx
                !y = (sin(x))**2
                y = (COMPLEX(x,-1.D0))**3
                array_x(i+1) = x
                array_y(i+1) = y
        end do
        !$omp end parallel do

        !r = REAL(trapezoid(array_x,array_y))
        r = REAL(trapz_parallel(array_x,array_y))
        write(6,'(2F12.6)') r, -5.D0/4.D0

        DEALLOCATE(array_x)
        DEALLOCATE(array_y)

END PROGRAM INTEGRATION_CHECKER

