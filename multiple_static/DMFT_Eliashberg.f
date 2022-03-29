      PROGRAM DMFT_Eliashberg
      INCLUDE 'mpif.h'
      !IMPLICIT NONE      
      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: smalleps = 1.e-5 
 
      INTEGER :: nq(3)  ! number of q points
      INTEGER :: iqx,iqy,iqz     ! index of q points
      INTEGER :: m      ! index of bosonic Matsubara frequency
      INTEGER :: n,nn   ! index of fermionic Matsubara frequency
      INTEGER :: n_shift, nn_shift, nn_shift_G

      REAL*8  :: T      ! temperature
      REAL*8  :: U      ! Hubbard U
      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency
      REAL*8  :: mu     ! chemical potential
      INTEGER :: nk(3), nkk(3)  ! number of k-point
      INTEGER :: numk, numkw

      REAL*8  :: k(3), kk(3)  ! k and k'
      INTEGER :: ikx, iky, ikz, ikkx, ikky, ikkz  ! index of k, k' points
      INTEGER :: ikwp_total
      INTEGER :: ik(3)

      REAL*8  :: mix       ! mixing parameter
      INTEGER :: num_iter  ! maximum number of iteration
      INTEGER :: iter      ! iteration number

      INTEGER :: r, nr                 ! number of "r" points
      INTEGER,ALLOCATABLE :: tran(:,:) ! R
      COMPLEX*16,ALLOCATABLE :: ham(:) ! H(R)
      COMPLEX*16,ALLOCATABLE :: sigi(:)! self-energy

      COMPLEX*16,ALLOCATABLE :: chi(:,:,:,:) ! irreducible susceptibility
      COMPLEX*16,ALLOCATABLE :: Green(:,:,:,:) ! Green function

      COMPLEX*16 :: Green_temp, Veff, chi_temp, chi_s, chi_c

      COMPLEX*16,ALLOCATABLE :: V_matrix(:,:), V_CT(:,:)
      COMPLEX*16,ALLOCATABLE :: tot_V_matrix(:,:)
      INTEGER*8 :: dim
      INTEGER :: r_index,c_index

      REAL*8,ALLOCATABLE :: ommesh(:)  ! fermionic Matsubara frequency

      REAL*8  :: initial ! initial guess
      LOGICAL :: negativem, negativenn
    

      REAL*8  :: numesh ! bosonic Matsubara frequency
      INTEGER :: i,i_shift,j  ! dummy index
      LOGICAL :: iffile

      ! for MPI
      INTEGER rank,total,ierr,errorcode
      INTEGER pr_proc,numkw_per,ikwp

      ! for lapack
      COMPLEX*16 :: leigenvalue      
      COMPLEX*16, ALLOCATABLE :: eigenvalues(:)
      REAL*8, ALLOCATABLE :: eigenvalues_real(:)
      INTEGER :: LDVL, LDVR, INFO, LWORK
      INTEGER :: location
      COMPLEX*16, ALLOCATABLE :: VL(:, :), VR(:, :), WORK(:)
      REAL*8, ALLOCATABLE :: RWORK(:) 

      WRITE(*,*) 'DMFT_Eliashberg starts:'

      INQUIRE(FILE = 'DMFT_Eliashberg.in', EXIST = iffile)
      IF (iffile .EQV. .false.) THEN
          WRITE(*,*) 'input is needed for DMFT_Eliashberg!'
          STOP
      ELSE
          OPEN(65, FILE = "DMFT_Eliashberg.in")

          READ(65,*) initial
          READ(65,*) U
          READ(65,*) T
          READ(65,*) mu
          READ(65,*) mix
          READ(65,*) num_iter

          READ(65,*) (nk(i),i=1,3)
          READ(65,*) Nc

          READ(65,*) nr
          ALLOCATE(tran(nr,3))
          READ(65,*) ((tran(i,j),j=1,3),i=1,nr)
          ALLOCATE(ham(nr))
          READ(65,*) (ham(i),i=1,nr)

          ALLOCATE(sigi(2*Nc))
          READ(65,*) (sigi(i),i=1,2*Nc)

          CLOSE(65)
      ENDIF !IF iffile

      ! here we 'nq' is equal 'nk'
      nq = nk
 
      ! here we 'nkk' is equal 'nk'    
      nkk = nk
 
      numk  = nk(1)*nk(2)*nk(3)
      numkw = nk(1)*nk(2)*nk(3)*(2*Nc)

      ALLOCATE(ommesh(2*Nc))
      ALLOCATE(chi(2*Nc, nq(1), nq(2), nq(3)))
      ALLOCATE(Green(Nc, nk(1), nk(2), nk(3)))

      dim = nk(1)*nk(2)*nk(3)*2*Nc
      ALLOCATE(V_matrix(nk(1)*nk(2)*nk(3)*2*Nc,nk(1)*nk(2)*nk(3)*2*Nc))
      ALLOCATE(V_CT(dim,dim))
      ALLOCATE(tot_V_matrix(dim,dim))

      LDVL = dim
      LDVR = dim
      ALLOCATE(eigenvalues(dim))
      ALLOCATE(eigenvalues_real(dim))
      ALLOCATE(VL(dim,dim))
      ALLOCATE(VR(dim,dim))

      LWORK = 2*dim
      ALLOCATE(WORK(LWORK))
      ALLOCATE(RWORK(LWORK))

      DO i = -Nc, Nc-1
          i_shift = i + (Nc+1)
          ommesh(i_shift) = (2*i+1)*pi*T
      ENDDO

      ! launch MPI
      CALL MPI_init(ierr)

      IF (ierr .NE. MPI_SUCCESS) THEN
          WRITE(*,*), 'Error starting MPI program. Terminating.'
          CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      ENDIF !IF ierr

      CALL MPI_Comm_size(MPI_Comm_World, total, ierr)
      CALL MPI_Comm_rank(MPI_Comm_World, rank, ierr)

      !here 'nq' is equal to 'nk' because we only calculate commensurate grids
      CALL do_chi(nk,T,Nc,mu,nk,ommesh,nr,tran,ham,sigi,chi,rank,total)
 
      CALL do_Green(nk,T,Nc,mu,ommesh,nr,tran,ham,sigi,Green,rank,total)

 
      ! pr_proc is the number of kpoints per CPU
      ! numk_per is equal to pr_proc except for the last CPU (highest
      ! rank)
      ! if nk/total is not an integer, numk_per = numk-rank*pr_proc

      pr_proc = FLOOR(numkw/DBLE(total)+0.999)


      numkw_per = pr_proc
      IF ((rank+1)*pr_proc .GT. numkw) THEN
          IF (numkw-rank*pr_proc .GT. 0) THEN
              numkw_per = numkw-rank*pr_proc
          ELSE
              numkw_per = 0
          ENDIF
      ENDIF

      WRITE(*,*) 'Rank: ', rank,  'pr_proc: ', pr_proc, ' numkw_per: ',
     1numkw_per

      !initialize V_matrix for all omega_n and k
      V_matrix = 0.0
      V_CT = 0.0

      Green_temp = 0.0

      WRITE(*,*) 'hello0'

      DO n = -Nc, Nc-1
          n_shift = n + (Nc+1)      
 
          DO ikx = 1, nk(1)
              DO iky = 1, nk(2)
                  DO ikz = 1, nk(3)

                      ! fill in the V_matrix
                      ! calculate the row index
                      r_index = 
     1                (n_shift-1)*nk(1)*nk(2)*nk(3) + 
     1                (ikx-1)*nk(2)*nk(3) + (iky-1)*nk(3) + ikz


                      ikwp = ikz+(iky-1)*nk(3)+(ikx-1)*nk(2)*nk(3)+
     1 (n_shift-1)*nk(1)*nk(2)*nk(3)-rank*pr_proc

                      ! parallel over different k points and omega
                      IF (ikwp .GT. 0 .AND. ikwp .LE. numkw_per) THEN

C                       IF (rank .eq. 0) THEN
C                           WRITE(*,*) 'ikwp: ', ikwp
C                       ENDIF
  
                          DO nn = -Nc, Nc-1

                              ! set 'negativem'
                              negativem = .FALSE.

                              ! get "m" and m_shift
                              m = n - nn

                              IF (m .lt. 0) THEN
                                  negativem = .TRUE.
                                  m = -m
                              ENDIF

                              m_shift = m + 1

                              ! set 'negativenn'

                              negativenn = .FALSE.

                              ! get nn_shift
                              nn_shift = nn + (Nc+1)
 
                              IF (nn .lt. 0) THEN
                                  negativenn = .TRUE.
                                  nn_shift_G = (ABS(nn)-1) + 1
                              ELSE
                                  nn_shift_G = nn + 1
                              ENDIF


                              DO ikkx = 1, nkk(1)
                                  DO ikky = 1, nkk(2)
                                      DO ikkz = 1, nkk(3)

                                          iqx = MOD(ikx - ikkx, nkk(1))
                                          iqy = MOD(iky - ikky, nkk(2))
                                          iqz = MOD(ikz - ikkz, nkk(3))

                                          IF (iqx .le. 0) iqx = 
     1                                    iqx + nkk(1)
                                          IF (iqy .le. 0) iqy = 
     1                                    iqy + nkk(2)
                                          IF (iqz .le. 0) iqz = 
     1                                    iqz + nkk(3)

                                          ! calculate chi_s and chi_c from chi

                                          IF (negativem .eqv. 
     1                                    .FALSE.) THEN
                                              chi_temp = 
     1                                        chi(m_shift,iqx,iqy,iqz)
                                          ELSE
                                              chi_temp = 
     1                                CONJG(chi(m_shift,iqx,iqy,iqz))
                                          ENDIF

                                          chi_s = 
     1                                    (chi_temp/(1.0-U*chi_temp))
                                          chi_c = 
     1                                    (chi_temp/(1.0+U*chi_temp))

                                          ! calculate Veff
                                          Veff = U+(3./2.)*U*U*chi_s-
     1                                    (1./2.)*U*U*chi_c

                                          IF (negativenn .eqv. 
     1                                    .FALSE.) THEN
                                              Green_temp = 
     1                                  Green(nn_shift_G,ikkx,ikky,ikkz)
                                          ELSE
                                              Green_temp =
     1                           CONJG(Green(nn_shift_G,ikkx,ikky,ikkz))
                                          ENDIF


C                                       WRITE(*,'(A, 6I3, 4F24.16)') 
C      1                                'YYY:', 
C      1                         nn, nn_shift, nn_shift_G, ikkx,ikky,ikkz,
C      1                         REAL(Green(nn_shift_G,ikkx,ikky,ikkz)),
C      1 AIMAG(Green(nn_shift_G,ikkx,ikky,ikkz)),REAL(Veff), AIMAG(Veff)

                                          ! fill in the V_matrix
C                                           r_index = 
C      1                                (n_shift-1)*nk(1)*nk(2)*nk(3) + 
C      1                                (ikx-1)*nk(2)*nk(3) + 
C      1                                (iky-1)*nk(3) + ikz

                                          c_index = 
     1                                (nn_shift-1)*nk(1)*nk(2)*nk(3) + 
     1                                (ikkx-1)*nk(2)*nk(3) + 
     1                                (ikky-1)*nk(3) + ikkz

                                          V_matrix(r_index,c_index) = 
     1                                Veff*(abs(Green_temp)**2)

C                                       WRITE(*,*) 'V(k,k_prime)',
C      1                                Veff*(abs(Green_temp)**2)

                                      ENDDO ! ikkz
                                  ENDDO ! ikky
                              ENDDO ! ikkx
                          ENDDO ! nn
                      ENDIF !ikwp
                  ENDDO ! ikz
              ENDDO ! iky
          ENDDO ! ikx
      ENDDO ! n

C       WRITE(*,*) Veff
C       WRITE(*,*) Green_temp

      V_matrix = -(T/numk)*V_matrix

C       WRITE(*,*) 'hello1'

      CALL MPI_REDUCE(V_matrix,tot_V_matrix,
     1dim*dim,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_Comm_World,ierr)

C       WRITE(*,*) 'hello2'
      
      IF (rank .eq. 0) THEN

          V_matrix = tot_V_matrix

          WRITE(*,*) 'V eff matrix'
          WRITE(*,*) V_matrix(1,1)
          WRITE(*,*) V_matrix(2,2)
          WRITE(*,*) V_matrix(3,3)
          WRITE(*,*) V_matrix(4,4)

          V_CT = CONJG(TRANSPOSE(V_matrix))

          V_CT = V_CT - V_matrix

          WRITE(*,*) 'V_matrix conjugate transpose'
          WRITE(*,*) 'imaginary part not 0'
C           DO i=1,dim
C               DO j=1,dim
C                   WRITE(*,*) i, j
C                   WRITE(*,*) V_CT(i,j)
C               ENDDO ! do j
C           ENDDO ! do i
          WRITE(*,*) 'difference between V and its conjugate transpose'
          WRITE(*,*) V_CT(1,1)
          WRITE(*,*) V_CT(2,2)
          WRITE(*,*) V_CT(3,3)
          WRITE(*,*) V_CT(4,4)

          CALL ZGEEV('N','N',dim,V_matrix,dim,eigenvalues,VL,LDVL,VR, 
     1    LDVR,WORK,LWORK,RWORK,INFO)

          eigenvalues_real = ABS(eigenvalues)
          location = MAXLOC(eigenvalues_real, 1)
          WRITE(*,*) 'location:', location
C           WRITE(*,*) 'eigenvalues:', eigenvalues

          leigenvalue = eigenvalues(location)

          OPEN(65, FILE = "DMFT_Eliashberg.out")             
          WRITE(65,'(4I3,2F24.16)') 2*Nc,nk(1),nk(2),nk(3),
     1    REAL(leigenvalue), AIMAG(leigenvalue)
C           DO i=1,dim
C               IF (AIMAG(eigenvalues(i)) .NE. 0) THEN
C                   WRITE(65,'(1I6,3F24.16)') i,REAL(eigenvalues(i)),
C      1            AIMAG(eigenvalues(i)),eigenvalues_real(i)
C               ENDIF
C           ENDDO
          CLOSE(65)

      ENDIF ! rank

      CALL MPI_Finalize(ierr)

      DEALLOCATE(chi)
      DEALLOCATE(Green)
      DEALLOCATE(ommesh)
      DEALLOCATE(V_matrix)
      DEALLOCATE(V_CT)
      DEALLOCATE(tot_V_matrix)
      DEALLOCATE(eigenvalues)
      DEALLOCATE(eigenvalues_real)
      DEALLOCATE(VL)
      DEALLOCATE(VR)
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)

      END PROGRAM DMFT_Eliashberg