PROGRAM DOS
        ! Calculates density of states using the tetrahedron interpolation.
        use declarations
        use compute
        use readio
        use OMP_LIB

        IMPLICIT NONE

        COMPLEX(DP), ALLOCATABLE        :: dist_evc_r_ik(:)
        INTEGER, ALLOCATABLE            :: array_tetra_index(:,:)
        REAL(DP), ALLOCATABLE           :: array_tetra(:,:),e_array(:,:)
        REAL(DP), DIMENSION(3,8)        :: c_m
        INTEGER, DIMENSION(8)           :: vertex
        COMPLEX(DP), ALLOCATABLE        :: func(:,:,:)
        REAL(DP), DIMENSION(3)          :: k, dk, test_k
        REAL(DP), DIMENSION(4)          :: array_energy, array_quantity
        REAL(DP), DIMENSION(3,4)        :: array_k
        INTEGER, ALLOCATABLE            :: index_array(:,:)

        CHARACTER(len=256)              :: fname

        REAL(DP)                        :: vol, sum_Fermi
        REAL(DP)                        :: c1, c2

        INTEGER                         :: iuwfcr,lrwfcr,ios,unf_recl,reclen
        INTEGER                         :: min_band, max_band, nbnd_active 
        INTEGER                         :: itetra,i,id,ik1,ik2,ik3,ik,iv,ibnd, loci
        LOGICAL                         :: itest

        REAL(DP)                        :: Efree, tolerance, dos_theory
        REAL(DP), ALLOCATABLE               :: DOS(:)
        INTEGER                         :: ie, nenergy

        itetra = 0
        nbnd_active = 0
        fname = trim(adjustl(trim(tmpdir) // trim(wfcrfile)))
        iuwfcr = 877
        ! lrwfc should have the same as the number written in the wfck2r output in the line 
        ! before the line "length of wfc in real space/per band = ..."
        ! The length of wfc in real space/per band should match lrwfcr*nks*2*8
        ! the factor 8 because DP precision
        ! factor 2 because complex number
        ! lrwfcr*nks*2*8*nbnd = total size of the binary file <prefix>.wfc_r
        ! that can be checked with ls -l
        lrwfcr = nr1*nr2*nr3  
        ios = 0
        vol = a1*a2*a3
        min_band = 1  ! intended to come down to the lowest active band
        max_band = 1     ! intended to go up to the highest active band
        dk(1) = 1.D0/n1   
        dk(2) = 1.D0/n2
        dk(3) = 1.D0/n3
        loci = 0
        sum_Fermi = 0.D0
        nenergy = 200
        ! set c1 and c2
        c1 = 0.5*sqrt(a1**2 + a2**2 + a3**2)
        c2 = 0.00
        tolerance = 1.D0/10000000.0

        ALLOCATE (dist_evc_r_ik(lrwfcr))
        ALLOCATE(array_tetra_index(8, 6*nbnd*nks))
        ALLOCATE(array_tetra(16, 6*nbnd*nks))
        ALLOCATE(index_array(2,nks))
        ALLOCATE(e_array(nbnd,nks))
        ALLOCATE(func(nr1,nr2,nr3))
        ALLOCATE(DOS(nenergy))
        
        
        !fill the e_array with free electron energy or tight binding energy
        do ik = 1, nks
            do ibnd = 1, 1
                k = indextok(ik, n1, n2, n3)
                !e_array(ibnd, ik) = sqrt((k(1)**2 +k(2)**2 + k(3)**2))
                e_array(ibnd, ik) = (k(1)**2 +k(2)**2 + k(3)**2)/2.0
                !e_array(ibnd, ik) = -2.D0*(cos(2*PI*k(1)) + cos(2*PI*k(2)) + cos(2*PI*k(3)))
                !write(6,'(3F8.3, 2X, F8.4)'), k, e_array(ibnd, ik)
            end do
        end do
       
        OPEN (unit = 91, file = "dos_free_parallel.dat", status="replace")
        OPEN (unit = 221, file = "tetra_list.dat", status="replace")
        
        func = 1.D0
        DOS = 0.D0

        DO ie = 1, nenergy 
            loci = 0
            array_tetra_index = (0, 0)
            array_tetra = (0.D0,0.D0)
            !Efree = -6.D0 + 12.D0*ie/float(nenergy)
            Efree = 0.0 + ie/float(2*nenergy)
            dos_theory = 2.D0*vol*2.D0**(3.D0/2)*sqrt(Efree)/(4.D0*PI**2)
            !dos_theory = 4.D0*PI*Efree**(2)
            itetra = 0
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
                                    id = ktoindex(k, n1,n2,n3)
                                    test_k = indextok(id, n1,n2,n3)

                                    if (SUM(ABS(k - test_k)) .GT. tolerance) then
                                        write(6,*) "index transformations are wrong"
                                         write(6,'(I5, 3F12.6, 2X, 3F12.6)') id,k,test_k
                                         write(6,*) "         "
                                        STOP 1
                                    end if
                                    !write(6, '(A, I5, 3F10.5)') "index and k-point:  ", &
                                    !        & ktoindex(k, n1, n2, n3), k
                                    ! create the corners of the submesh cell starting from k. 
                                    CALL scan_BZ(k, dk, c_m)
                                    itetra = itetra + 6
                                    ! One consistency check 
                                    if ( (ktoindex(c_m(:,6),n1,n2,n3) < 1) .OR. &
                                            & (ktoindex(c_m(:,6), n1,n2,n3) > nks) ) then
                                            write(6, '(A)') " Error in sampling the BZ. k reaches beyond the grid"
                                            write(6,'(A, I8, 3F8.4, I8, 3F8.4)')"corners:indices and k", &
                                                & ktoindex(c_m(:,3),n1,n2,n3), c_m(:,3), &
                                                & ktoindex(c_m(:,6),n1,n2,n3), c_m(:,6)
                                            STOP 1
                                    end if
                                    ! check if any of the 6 tetrahedrons, for any of the bands, cross the Fermi level. 
                                    !If it does, ! fill the array_tetra and array_tetra_index up
                                    call chk_tetra(c_m(:,1),c_m(:,2),c_m(:,3),c_m(:,6),&
                                    & e_array,min_band,max_band,array_tetra_index,array_tetra,Efree,loci)
                                    call chk_tetra(c_m(:,2),c_m(:,3),c_m(:,4),c_m(:,6),&
                                    & e_array,min_band,max_band,array_tetra_index,array_tetra,Efree,loci)
                                    call chk_tetra(c_m(:,1),c_m(:,3),c_m(:,5),c_m(:,6),&
                                    & e_array,min_band,max_band,array_tetra_index,array_tetra,Efree,loci)
                                    call chk_tetra(c_m(:,3),c_m(:,4),c_m(:,6),c_m(:,8),&
                                    & e_array,min_band,max_band,array_tetra_index,array_tetra,Efree,loci)
                                    call chk_tetra(c_m(:,3),c_m(:,5),c_m(:,6),c_m(:,7),&
                                    & e_array,min_band,max_band,array_tetra_index,array_tetra,Efree,loci)
                                    call chk_tetra(c_m(:,3),c_m(:,6),c_m(:,7),c_m(:,8),&
                                    & e_array,min_band,max_band,array_tetra_index,array_tetra,Efree,loci)
    
                            end do
                    end do
            end do
            write(221,'(A,F5.3,A,I5,A,I5)') "Energy ", Efree, & 
                    & " active tetrahedrons ", loci, " total tetrahedrons ", itetra
            ! loop over as many active tetrahedrons as you have
            DO itetra = 1, loci
                ! get the energies of the 4 corners
                array_energy = array_tetra(13:16, itetra)
                ! Loop over the vertices of each tetrahedron
                DO iv = 0,3
                    ! get the index of the k-vector for a single vertex
                    ik = array_tetra_index(iv+1, itetra)
                    ! Get the corresponding band
                    ibnd = array_tetra_index(iv+5, itetra)
                    ! get the corresponding k-vector (can get it from indextok as well)
                    array_k(:,iv+1) = array_tetra(3*iv+1:3*iv+3, itetra)
                    !array_quantity(iv + 1) =  f_kkp(ik, ibnd, ik, ibnd)
                    array_quantity(iv+1) = 1.D0
                END DO
                DOS(ie) = DOS(ie) + tetra(array_k, array_energy, array_quantity, Efree) 
            END DO
            !write(91, '(F8.4, 2X, F8.4, 2X, F8.4)') Efree, DOS(ie), dos_theory
        END DO


        DO ie = 1, nenergy
            Efree = 0.0 + ie/float(2*nenergy)
            dos_theory = 2.D0*vol*2.D0**(3.D0/2)*sqrt(Efree)/(4.D0*PI**2)
            write(91, '(F8.4, 2X, F8.4, 2X, F8.4)') Efree, DOS(ie), dos_theory
        END DO     

        DEALLOCATE (dist_evc_r_ik)
        DEALLOCATE(array_tetra_index)
        DEALLOCATE(array_tetra)
        DEALLOCATE(func)
        DEALLOCATE(index_array)
        DEALLOCATE(e_array)
        DEALLOCATE(DOS)
        CLOSE(unit=221)
        CLOSE(unit=91)
        CLOSE(iuwfcr) 

END PROGRAM DOS

