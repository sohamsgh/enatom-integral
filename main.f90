PROGRAM MAIN
        ! Calculates the Fermi surface integral
        use declarations
        use compute
        use readio

        IMPLICIT NONE

        COMPLEX(DP), ALLOCATABLE        :: dist_evc_r_ik(:), dist_evc_r_ikp(:)
        INTEGER, ALLOCATABLE            :: array_tetra_index(:,:)
        REAL(DP), ALLOCATABLE           :: array_tetra(:,:),e_array(:,:)
        REAL(DP), DIMENSION(3,8)        :: c_m
        INTEGER, DIMENSION(8)           :: vertex
        REAL(DP), ALLOCATABLE           :: iearray(:)
        REAL(DP), ALLOCATABLE           :: f_kkp(:,:,:,:)
        COMPLEX(DP), ALLOCATABLE        :: func(:,:,:)
        REAL(DP), DIMENSION(3)          :: k, dk
        REAL(DP), DIMENSION(4)          :: array_energy,array_energyp,array_quantity
        REAL(DP), DIMENSION(4)          :: g_k
        REAL(DP), DIMENSION(3,4)        :: array_k,array_kp
        INTEGER, ALLOCATABLE            :: index_array(:,:)

        CHARACTER(len=256)              :: fname

        REAL(DP)                        :: vol, sum_Fermi
        REAL(DP)                        :: c1, c2

        INTEGER                         :: nrec_ik,nrec_ikp
        INTEGER                         :: iuwfcr,lrwfcr,ios,unf_recl,reclen,ios_ik,ios_ikp 
        INTEGER                         :: min_band, max_band, nbnd_active 
        INTEGER                         :: itetra,itetrap,ik1,ik2,ik3,ik,ikp,iv,ivp,ibnd,ibndp, loci, ic2 
        LOGICAL                         :: itest

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
        min_band = nbnd  ! intended to come down to the lowest active band
        max_band = 1     ! intended to go up to the highest active band
        dk(1) = 1.0/n1   
        dk(2) = 1.0/n2
        dk(3) = 1.0/n3
        loci = 0
        sum_Fermi = 0.0

        ! set c1 and c2
        c1 = 0.5*sqrt(a1**2 + a2**2 + a3**2)
        c2 = 0.00
        !call READINPUT(c1, c2)
        !write(6,*) c1, c2

        ALLOCATE (dist_evc_r_ik(lrwfcr))
        ALLOCATE (dist_evc_r_ikp(lrwfcr))
        ALLOCATE(array_tetra_index(8, 6*nbnd*nks))
        ALLOCATE(array_tetra(16, 6*nbnd*nks))
        ALLOCATE(index_array(2,nks))
        ALLOCATE(e_array(nbnd,nks))
        ALLOCATE(iearray(nks))
        ALLOCATE(f_kkp(nks,nbnd,nks,nbnd))
        ALLOCATE(func(nr1,nr2,nr3))
        array_tetra_index = (0, 0)
        array_tetra = (0.D0,0.D0)

        unf_recl = 2*DIRECT_IO_FACTOR*lrwfcr
        inquire(iolength=reclen) dist_evc_r_ik
        if (reclen .NE. unf_recl) then
                write(6, *) "Record length in binary wfc file is not correct"
                STOP 1
        end if
        
        ! Run the READEVAL routine to get the eigenvalues for each k-point and band
        ! We don't need index_array for the tetrahedron method 
        CALL READEVAL(del, index_array, e_array)

        ! Get the lowest band and highest bands that cross the Fermi level
        ! and look at only those bands.
        ! It is possible that one band entirely below the fermi level can connect to 
        ! Another band entirely above the Fermi surface by being two corners of a tetrahedron. 
        ! But those should be rare and the contribution of such bands to the final integral 
        ! should be small enough to ignore. 
        do ibnd = 1, nbnd
                ! Get the energy for each band for all k-points
                iearray = e_array(ibnd,:)
                if ( (minval(iearray) .LT. EFermi) .AND. (maxval(iearray) .GT. EFermi) ) then
                        nbnd_active = nbnd_active + 1
                        ! #min_band comes down from nbnd to the first match and stays there.
                        if (ibnd .LE. min_band) then 
                                min_band = ibnd
                        end if
                        ! #max_band starts at 1 and goes up till the outer loop in entered.
                        if (ibnd .GE. max_band) then
                                max_band = ibnd
                        end if
                end if
        end do
        if (nbnd_active == 0) then
                write(6,*) "No band crosses Fermi surface, &
                        &can't do a Fermi surface integral"
        else
                !write(6, *) "min_band:", min_band
                !write(6, *) "max_band:", max_band
        end if
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
                                !write(6, '(A, I5, 3F10.5)') "index and k-point:", &
                                        !& ktoindex(k, n1, n2, n3), k
                                ! create the corners of the submesh cell starting from k. 
                                CALL scan_BZ(ik1, ik2, ik3, k, dk, c_m)
                                do ik = 1,8
                                        vertex(ik) = ktoindex(c_m(:,ik), n1, n2, n3)
                                end do
                                itetra = itetra + 1
                                ! Write the vertices
                                !write(6,'(A,I6)') "tetra", itetra
                                !write(6,'(I5)') vertex(1)
                                !write(6,'(I5)') vertex(2)
                                !write(6,'(I5)') vertex(3)
                                !write(6,'(I5)') vertex(6)
                                itetra = itetra + 1
                                !write(6,'(A,I6)') "tetra", itetra
                                !write(6,'(I5)') vertex(2)
                                !write(6,'(I5)') vertex(3)
                                !write(6,'(I5)') vertex(4)
                                !write(6,'(I5)') vertex(6)
                                itetra = itetra + 1
                                !write(6,'(A,I6)') "tetra", itetra
                                !write(6,'(I5)') vertex(1)
                                !write(6,'(I5)') vertex(3)
                                !write(6,'(I5)') vertex(5)
                                !write(6,'(I5)') vertex(6)
                                itetra = itetra + 1
                                !write(6,'(A,I6)') "tetra", itetra
                                !write(6,'(I5)') vertex(3)
                                !write(6,'(I5)') vertex(4)
                                !write(6,'(I5)') vertex(6)
                                !write(6,'(I5)') vertex(8)
                                itetra = itetra + 1
                                !write(6,'(A,I6)') "tetra", itetra
                                !write(6,'(I5)') vertex(3)
                                !write(6,'(I5)') vertex(5)
                                !write(6,'(I5)') vertex(6)
                                !write(6,'(I5)') vertex(7)
                                itetra = itetra + 1
                                !write(6,'(A,I6)') "tetra",itetra
                                !write(6,'(I5)') vertex(3)
                                !write(6,'(I5)') vertex(6)
                                !write(6,'(I5)') vertex(7)
                                !write(6,'(I5)') vertex(8)


                                ! One consistency check 
                                if ( (ktoindex(c_m(:,6),n1,n2,n3) < 1) .OR. &
                                        & (ktoindex(c_m(:,6), n1,n2,n3) > nks) ) then
                                        write(6, '(A, 3F6.3)') " Error in sampling the BZ. k &
                                                & reaches beyond the grid at", c_m(:,6)
                                        STOP 1
                                end if
                                ! check if any of the 6 tetrahedrons, for any of the bands, cross the Fermi level. 
                                !If it does, ! fill the array_tetra and array_tetra_index up
                                call chk_tetra(c_m(:,1),c_m(:,2),c_m(:,3),c_m(:,6),&
                                & e_array,min_band,max_band,array_tetra_index,array_tetra,EFermi,loci)
                                call chk_tetra(c_m(:,2),c_m(:,3),c_m(:,4),c_m(:,6),&
                                & e_array,min_band,max_band,array_tetra_index,array_tetra,EFermi,loci)
                                call chk_tetra(c_m(:,1),c_m(:,3),c_m(:,5),c_m(:,6),&
                                & e_array,min_band,max_band,array_tetra_index,array_tetra,EFermi,loci)
                                call chk_tetra(c_m(:,3),c_m(:,4),c_m(:,6),c_m(:,8),&
                                & e_array,min_band,max_band,array_tetra_index,array_tetra,EFermi,loci)
                                call chk_tetra(c_m(:,3),c_m(:,5),c_m(:,6),c_m(:,7),&
                                & e_array,min_band,max_band,array_tetra_index,array_tetra,EFermi,loci)
                                call chk_tetra(c_m(:,3),c_m(:,6),c_m(:,7),c_m(:,8),&
                                & e_array,min_band,max_band,array_tetra_index,array_tetra,EFermi,loci)

                        end do
                end do
        end do
        !write(6,*) "loci =", loci
        OPEN (unit = 221, file = "FSfile.dat", status="replace")
        write(221,'(A)') "Total number ofo tetrahedra"
        write(221,'(A, I5)') "Total number of tetrahedrons", loci
        do ik = 1, loci
                write(221,'(A, I5)') "Tetrahedrons", ik
                write(221, *) array_tetra(1:3,ik)
                write(221, *) array_tetra(4:6,ik)
                write(221, *) array_tetra(7:9,ik)
                write(221, *) array_tetra(10:12,ik)
                write(221,'(A)') "        "
        end do
        CLOSE(unit=221)

        ! open binary wavefunction file
        ! First lets test the tetrahedron interpolation routine. For that, 
        ! Set the input quantity to 1.0.

        ! reads the real space wavefunctions from the 
        ! <prefix>.wfc_r file. This is a binary file
        ! lrwfcr = nr1*nr2*nr2 is the record length and nrec is the record location.

        OPEN (unit = iuwfcr, file = fname, iostat = ios, form = 'unformatted', &
                status = 'old', action='read', access = 'direct', recl = reclen)
        ! Calculate the weight of the double dirac delta , which is
        ! |<psi_(k', n') (x,y,z))| f(x,y,z)| psi_(k,n)(x,y,z)>|^2 
        ! for all k, n, k', n'
        ! This is faster than doing them inside the loops of the tetrahedron 
        ! Because in the tetrahedron list there are repetitions of k and band.

        !    ############################################################################
        !
        !                               Some testing
        !
        !    ############################################################################
        
        ! The first job is to test the tetrahedron formulation
        ! One test is to calculate DOS(E) for a free electron gas
        itest = .FALSE. 
        if ( itest .EQV. .TRUE.) then
            do energy = 1, 
            do ik = 1, nks
                do ibnd = min_band, max_band
                    nrec_ik = (ik-1)*nbnd + ibnd
                    dist_evc_r_ik = (0.0, 0.0)

                end do
            end do
        end if




      do ic2 = 1, 30
        c2 = 0.01 + ic2*(3.0-0.1)/30

        func = 14*READENTATOM(c1,c2)
        do ik = 1, nks
            do ikp = 1, nks
                do ibnd = min_band, max_band
                    do ibndp = min_band, max_band
                        nrec_ik = (ik-1)*nbnd+ibnd
                        nrec_ikp = (ikp-1)*nbnd+ibndp
                        dist_evc_r_ik = (0.0,0.0)
                        dist_evc_r_ikp = (0.0,0.0)
                        ! Read the wavefunction for a given k and n
                        ! Its a 1-D array, the flatted form of psi(x, y, z)
                        ! Assume x moves the fastest and z the slowest
                        READ( unit = iuwfcr, REC = nrec_ik, IOSTAT = ios_ik ) dist_evc_r_ik
                        ! Read the wavefunction for a given k and n
                        READ( unit = iuwfcr, REC = nrec_ikp, IOSTAT = ios_ikp ) dist_evc_r_ikp
                        ! Calculate the square of the integral 
                        f_kkp(ik,ibnd,ikp,ibndp) = expect(dist_evc_r_ik, dist_evc_r_ikp, func)
                    end do
                end do
            end do
        end do    


        sum_Fermi = 0.0
        ! loop over as many active tetrahedrons as you have
        DO itetra = 1, loci
            g_k = (0.0)
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
                ! loop over as many active tetrahedrons as you have
                DO itetrap = 1, loci
                    array_energyp = array_tetra(13:16, itetrap)
                    ! Loop over the vertices of each tetrahedron
                    DO ivp = 0,3
                        ikp = array_tetra_index(iv+1, itetra)   
                        ibndp = array_tetra_index(iv+5, itetra)   
                        array_kp(:,ivp+1) = array_tetra(3*ivp+1:3*ivp+3, itetrap)
                        array_quantity(ivp + 1) =  f_kkp(ik, ibnd, ikp, ibndp)
                    END DO    
                    g_k(iv + 1) = g_k(iv + 1) + tetra(array_kp, array_energyp, array_quantity, EFermi)
                END DO
            END DO
            sum_Fermi = sum_Fermi + tetra(array_k, array_energy, g_k, EFermi) 
        END DO
        sum_Fermi = sum_Fermi/(vol*nks)
        write(6,*) c2, sum_Fermi
        end do

        DEALLOCATE (dist_evc_r_ik)
        DEALLOCATE (dist_evc_r_ikp)
        DEALLOCATE(f_kkp)
        DEALLOCATE(array_tetra_index)
        DEALLOCATE(array_tetra)
        DEALLOCATE(iearray)
        DEALLOCATE(func)
        DEALLOCATE(index_array)
        DEALLOCATE(e_array)


END PROGRAM MAIN
