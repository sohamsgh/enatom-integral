PROGRAM INDEX_CHECKER
        ! Calculates the Fermi surface integral
        use declarations
        use compute

        IMPLICIT NONE

        REAL(DP), DIMENSION(3)          :: k
        INTEGER                         :: id

        write(6,'(A)') "# ktoindex(k, n1, n2, n3) result(id)"
        write(6,'(A)') "#     id   Kx      ky      kz"
        do id = 1, nks
                k = indextok(id, n1, n2, n3)
                if (id .NE. ktoindex(k, n1, n2, n3)) then
                        write(6,*) "Houston we have a problem"
                else
                        write(6,'(I8, 3F9.4)') id, k
                end if
        end do

END PROGRAM INDEX_CHECKER

