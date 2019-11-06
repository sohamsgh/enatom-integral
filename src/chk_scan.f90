PROGRAM chk_scan
        ! checks the generation of tetrahedrons by scanning the BZ.
        use declarations
        use compute
        use readio

        IMPLICIT NONE

        REAL(DP), DIMENSION(3,8)        :: c_m1, c_m2, c_m_diff
        REAL(DP), DIMENSION(3,4)        :: vertices
        REAL(DP), DIMENSION(3)          :: k, dk, a, b, c, d, cp, n

        INTEGER                         :: ik1,ik2,ik3,ik,iv,iloc
        REAL(DP)                        :: tot1, tot2, tetra_vol, known_vol, tolerance

        dk(1) = 1.0/n1   
        dk(2) = 1.0/n2
        dk(3) = 1.0/n3
        tot1 = 0.0
        tot2 = 0.0
        a = 0.0
        b=0.0
        c=0.0
        d=0.0
        n(1) = n1
        n(2) = n2
        n(3) = n3
        tolerance = 1.D0/1000000000
        ! Now the scanning of the BZ to create tetrahedrons
        ! Scan from the extreme negative corner to the extreme positive corner 
        ! Stop at the penultimate grid point. 
        ! We will treat the boundary points seperately. 
        !write(6, '(A)') "Building tetrahedrons"
        !write(6, '(A)') "----------------------------"
        do ik1 = -(n1/2), (n1-1)/2
                do ik2 = -(n2/2), (n2-1)/2
                        do ik3 = -(n3/2), (n3-1)/2
                                k = (/ REAL(ik1)/n1, REAL(ik2)/n2, REAL(ik3)/n3 /)
                                CALL scan_BZ(k, dk, c_m1)
                                CALL scan_BZ2(ik1, ik2, ik3, k, dk, c_m2)
                                ! One consistency check 
                                if ( (ktoindex(c_m1(:,6),n1,n2,n3) < 1) .OR. &
                                        & (ktoindex(c_m1(:,6), n1,n2,n3) > nks) ) then
                                        write(6, '(A, 3F6.3)') " Error in sampling the BZ. k &
                                                & reaches beyond the grid at", c_m1(:,6)
                                        STOP 1
                                end if
                                ! Check if scan_BZ and scan_BZ2 both does the same thing.
                                c_m_diff = ABS(c_m1 - c_m2)
                                do iv = 1, 8
                                        tot1 = tot1 + sum(c_m_diff(:,iv))
                                end do
                                tot2 = tot2 + tot1 
                                !write(6,'(A, F12.6)') "Should be zero ", tot1
                                ! Check volume of tetrahedron based on volume
                                !https://en.wikipedia.org/wiki/Tetrahedron#Volume
                                d = c_m2(:,3)
                                a = c_m2(:,1) - d
                                b = c_m2(:,2) - d
                                c = c_m2(:,6) - d
                                d = 0.0
                                vertices(:,1) = d
                                vertices(:,2) = a
                                vertices(:,3) = b
                                vertices(:,4) = c
                                do iv = 1, 4
                                        do iloc = 1,3
                                                if (vertices(iloc, iv) .GT. REAL(n(iloc)-2)/n(iloc)) then
                                                        vertices(iloc,iv) = vertices(iloc,iv) - &
                                                                & REAL(n(iloc) -2)/n(iloc)
                                                else if (vertices(iloc, iv) .LT. -REAL(n(iloc)-2)/n(iloc)) then
                                                        vertices(iloc,iv) = vertices(iloc,iv) + &
                                                                & REAL(n(iloc) -2)/n(iloc)
                                                else
                                                end if
                                        end do
                                end do
                                cp = cross(vertices(:,3), vertices(:,4))
                                tetra_vol = ABS(vertices(1,2)*cp(1) + vertices(2,2)*cp(2) + vertices(3,2)*cp(3))/6.D0
                                if ((tetra_vol - ABS(dk(1)*dk(2)*dk(3))/6) .GT. tolerance) then 
                                        write(6,'(A, F12.8)') "Should be zero", &
                                                & tetra_vol, ABS(dk(1)*dk(2)*dk(3))/6.D0
                                end if
                        end do
                end do
        end do
        write(6,'(A, F12.6)') "Total Should be zero ", tot2
        

END PROGRAM chk_scan
