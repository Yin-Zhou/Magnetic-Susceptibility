      PROGRAM do_EP_super
      INCLUDE 'mpif.h'
      !IMPLICIT NONE
      !This program calculate the spectral function of a non-interacting
      !system in the superconducting state with electron-phonon coupling.

      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: smalleps = 1.e-5

      ! declare variables
      ! ###############################################################
      
      COMPLEX*16,DIMENSION(2,2) :: tau0,tau1,tau2,tau3      

      INTEGER :: numq  ! number of q points
      INTEGER :: numk ! number of k-point
      INTEGER :: numnu
      INTEGER :: iqx,iqy,iqz     ! index of q points
      INTEGER :: norb     ! number of orbital (currently it is 1)
      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency
      INTEGER :: Ncq, inu

      REAL*8 :: T !temperature
      REAL*8 :: eps_acustic !cutoff for acoustic phonons
      REAL*8 :: strength    !an artificial strength 
      REAL*8 :: deltakq 

      REAL*8,ALLOCATABLE :: ommesh(:)  ! fermionic Matsubara frequency
      REAL*8,ALLOCATABLE :: ommeshq(:)
      REAL*8,ALLOCATABLE :: kpoints(:,:),qpoints(:,:)

      COMPLEX*16,ALLOCATABLE :: Gv(:,:,:)
      COMPLEX*16,ALLOCATABLE :: temp1(:,:)
      COMPLEX*16 :: temp2,temp3
      COMPLEX*16,ALLOCATABLE :: temp4(:,:),temp5(:,:),tot_temp5(:,:)
      COMPLEX*16,ALLOCATABLE :: GREEN(:,:,:,:)
      REAL*8,ALLOCATABLE :: spectral(:,:)

      INTEGER :: ikx, iky, ikz  ! index of k, k' points

      REAL*8,ALLOCATABLE :: phonon_frequency(:,:)
      REAL*8,ALLOCATABLE :: enk(:),enkq(:,:)

      INTEGER :: i,j,l,m,n ! dummy index

      ! for MPI
      INTEGER rank,total,ierr,errorcode
      INTEGER pr_proc,numkw_per,ikwp
      ! for LAPACK
      INTEGER :: info
      INTEGER,DIMENSION(2) :: ipiv
      COMPLEX*16,DIMENSION(2) :: work

      ! declare variables finished
      ! ###############################################################

      ! get input
      ! ###############################################################
      WRITE(*,*) 'do_EP_super starts:'

      INQUIRE(FILE = 'EP.in', EXIST = iffile)
      IF (iffile .EQV. .false.) THEN
          WRITE(*,*) 'input is needed for do_EP_super!'
          STOP
      ELSE
          OPEN(65, FILE = "EP.in")

          READ(65,*) T

          READ(65,*) norb

          READ(65,*) numk
          READ(65,*) Nc

          READ(65,*) numq
          READ(65,*) Ncq

          READ(65,*) numnu

          !read in enk, enk+q
          ALLOCATE(enk(numk))
          READ(65,*) (enk(i),i=1,numk)

          ALLOCATE(enkq(numk,numq))
          READ(65,*) ((enkq(i,j),j=1,numq),i=1,numk)

          !read in k-points (along a path)
          ALLOCATE(kpoints(numk,3))
          READ(65,*) ((kpoints(i,j),j=1,3), i=1,numk)

          !read in q-points (high-symmetry qpoints)
          ALLOCATE(qpoints(numq,3))
          READ(65,*) ((qpoints(i,j),j=1,3), i=1,numq)

          !read in phonon frequency w_qv
          ALLOCATE(phonon_frequency(numq,numnu))
          READ(65,*) ((phonon_frequency(i,j),j=1,numnu),i=1,numq)

          !read in electron-phonon matrix g_v(k,q)
          ALLOCATE(Gv(numk,numq,numnv))
          READ(65,*) (((Gv(i,j,l),l=1,nu),j=1,numq),i=1,numk)       

          CLOSE(65)
      ENDIF !IF iffile

      ! get input finished
      ! ###############################################################
      
      ! allocate variables
      ! ###############################################################
      ALLOCATE(ommesh(2*Nc))
      ALLOCATE(ommeshq(2*Ncq))
      ALLOCATE(GREEN(numk,2*Nc,2,2))
      ALLOCATE(spectral(numk,2*Nc))

      ALLOCATE(temp1(2,2))
      ALLOCATE(temp4(2,2))
      ALLOCATE(temp5(2,2))
      ALLOCATE(tot_temp5(2,2))
      ! allocate variables finished
      ! ###############################################################

      ! introduce Pauli matrices

      tau0 = reshape((/1.0,0.0,0.0,1.0/),shape(tau0))
      tau1 = reshape((/0.0,1.0,1.0,0.0/),shape(tau1))
      tau2 = reshape((/(0.0,0.0),(0.0,1.0),(0.0,-1.0),(0.0,0.0)/),shape(tau2))
      tau3 = reshape((/1.0,0.0,0.0,-1.0/),shape(tau3))

      DO i = -Nc, Nc-1
          i_shift = i + (Nc+1)
          ommesh(i_shift) = (2*i+1)*pi*T
      ENDDO

      DO i = -Ncq, Ncq-1
          i_shift = i + (Ncq+1)
          ommeshq(i_shift) = 2*i*pi*T
      ENDDO

      ! launch MPI
      CALL MPI_init(ierr)

      IF (ierr .NE. MPI_SUCCESS) THEN
          WRITE(*,*), 'Error starting MPI program. Terminating.'
          CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      ENDIF !IF ierr

      CALL MPI_Comm_size(MPI_Comm_World, total, ierr)
      CALL MPI_Comm_rank(MPI_Comm_World, rank, ierr)
 
      ! pr_proc is the number of kpoints per CPU
      ! numk_per is equal to pr_proc except for the last CPU (highest
      ! rank)
      ! if nk/total is not an integer, numk_per = numk-rank*pr_proc

      pr_proc = FLOOR(numq/DBLE(total)+0.999)

      numq_per = pr_proc
      IF ((rank+1)*pr_proc .GT. numq) THEN
          IF (numq-rank*pr_proc .GT. 0) THEN
              numq_per = numq-rank*pr_proc
          ELSE
              numq_per = 0
          ENDIF
      ENDIF

      sigi = 0.0
      GREEN = 0.0

      WRITE(*,*) 'Rank: ', rank,  'pr_proc: ', pr_proc, ' numkw_per: ', numkw_per
      
      ! loop over k frequency
      DO n = -Nc, Nc-1
          n_shift = n + (Nc+1)
          ! loop over k points
          DO ikx = 1, numk

              ! G = i w tau0 - ep tau3 - delta tau1
              ! delta0 = 35 meV
              temp1=(0.0,1.0)*ommesh(n_shift)*tau0-enk(ikx)*tau3-(35*(COS(kpoints(ikx,1))-COS(kpoints(ikx,2)))/2)*tau1

              temp5 = 0.0
              tot_temp5 = 0.0
              ! loop over high-symmetry qpoints
              DO iqx = 1, numq

                  ikwp = iqx-rank*pr_proc
                  !parallel over different q points 
                  IF (ikwp .GT. 0 .AND. ikwp .LE. numq_per) THEN    
                          
                      ! loop over q frequency
                      DO m = -Ncq, Ncq-1
                          m_shift = m + (Ncq+1)
                          ! loop over nu
                          DO inu = 1,numnu

                              temp2=Gv(ikx,iqx,inu)**2

                              temp3=2*phonon_frequency(iqx,inu)/(((0.0,1.0)*ommeshq(m_shift))**2-phonon_frequency(iqx,inu)**2)

                              temp4=((0.0,1.0)*ommesh(n_shift)+(0.0,1.0)*ommeshq(m_shift))*tau0

                              deltakq=35*(COS(kpoints(ikx,1)+qpoints(iqx,1))-COS(kpoints(ikx,2)+qpoints(iqx,2)))/2

                              temp4=temp4+enkq(ikx,iqx)*tau3-deltakq*tau1

                              temp4=temp4/(((0.0,1.0)*ommesh(n_shift)+(0.0,1.0)*ommeshq(m_shift))**2-enkq(ikx,iqx)**2-deltakq**2)

                              temp5 = temp5 + temp2*temp3*temp4
                          ENDDO ! inu
                      ENDDO ! m
                  ENDIF !ikwp       
              ENDDO ! iqx 
              ! use MPI to do the sum over q
              CALL MPI_REDUCE(temp5,tot_temp5,4,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_Comm_World,ierr)

              IF (rank .EQ. 0) THEN
                  GREEN(ikx,n_shift,:,:)=temp1-T*tot_temp5
                  CALL ZGETRF(2,2,GREEN(ikx,n_shift,:,:),2,ipiv,info)
                  CALL ZGETRI(2,GREEN(ikx,n_shift,:,:),2,ipiv,work,2,info)
                  spectral(ikx,n_shift)=-AIMAG(GREEN(ikx,n_shift,1,1))/pi
              ENDIF ! rank
          ENDDO ! ikx
      ENDDO ! n

      ! output
      IF (rank .EQ. 0) THEN
          OPEN(66, FILE = "EP.out")             
      
          WRITE(66, '(A,I4)') "number of k points: ", numk
          WRITE(66, '(A,I4)') "number of omega: ", Nc
          WRITE(66, '(A)') "K-points are: "
          DO i = 1, numk
              WRITE(66,'(3F24.16)') kpoints(i,1), kpoints(i,2), kpoints(i,3)
          ENDDO

          CLOSE(66)

          !write down spectral function

          distance = 0.0

          OPEN(67, FILE= "spectral.out")
  
              DO ikx = 1, numk
                  DO n = 1, 2*Nc
                      WRITE(67,'(2F10.6, F24.16)') distance, ommesh(n), spectral(ikx,n)
                  ENDDO ! n               
                  IF (iq .lt. numq) THEN
                      distance = distance + sqrt((kpoints(ikx+1,1) - kpoints(ikx,1))**2 +(kpoints(ikx+1,2) - kpoints(ikx,2))**2 + (kpoints(ikx+1,3) - kpoints(ikx,3))**2)
                  ENDIF
              ENDDO !numk

          CLOSE(67)

      ENDIF ! rank      

      CALL MPI_Finalize(ierr)      

      DEALLOCATE(enk)
      DEALLOCATE(enkq)
      DEALLOCATE(kpoints)
      DEALLOCATE(qpoints)
      DEALLOCATE(phonon_frequency)
      DEALLOCATE(Gv)
      DEALLOCATE(ommesh)
      DEALLOCATE(ommeshq)
      DEALLOCATE(GREEN)
      DEALLOCATE(spectral)
      DEALLOCATE(temp1)
      DEALLOCATE(temp4)
      DEALLOCATE(temp5)
      DEALLOCATE(tot_temp5)

      END PROGRAM do_EP_super