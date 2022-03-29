      PROGRAM do_EP_self
      INCLUDE 'mpif.h'
      !IMPLICIT NONE
      !This program calculate the band-resolved self-energy of a non-interacting
      !system in the normal state arising from electron-phonon interaction.
     
 
      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: eps= 1.e-7

      ! declare variables
      ! ###############################################################

      INTEGER :: numomega ! number of  omega (omega is real)
      INTEGER :: numk     ! number of k-points
      INTEGER :: numq     ! number of q-points
      INTEGER :: numnv    ! number of phonon branch
      INTEGER :: norb     ! number of orbital (currently it is 1)     
      INTEGER :: numr, r  ! number of R points and dummy index for R

      !INTEGER,ALLOCATABLE :: tran(:,:) ! for R vectors
      REAL*8 :: epsk, epskq !epsilon_{k} and {k+q} 
      REAL*8 :: wqnv        !omega_{qv}
      REAL*8 :: mu          !chemical potential
      REAL*8 :: T           !temperature
      REAL*8 :: delta       !smearing width
      REAL*8 :: eps_acustic !cutoff for acoustic phonons

      !COMPLEX*16,ALLOCATABLE :: ham(:)  !H(R)

      REAL*8,ALLOCATABLE :: kpoints(:,:)   ! k-path
      REAL*8,ALLOCATABLE :: qpoints(:,:)   ! high-symmetry q points

      COMPLEX*16,ALLOCATABLE :: Gv(:,:,:)  ! electron-phonon matrix
      COMPLEX*16,ALLOCATABLE :: phonon_frequency(:,:) ! w_qv
      COMPLEX*16,ALLOCATABLE :: sigi(:)  ! sigi depends on k  
      COMPLEX*16,ALLOCATABLE :: tot_sigi(:)
      COMPLEX*16 :: temp1,temp2,temp3
      COMPLEX*16,ALLOCATABLE :: Green(:,:)
      REAL*8,ALLOCATABLE :: spectral(:,:)
      !REAL*8,ALLOCATABLE :: tot_spectral(:,:)
      REAL*8 :: nb,nf    ! FD and BE occupancy
      REAL*8,ALLOCATABLE :: weights(:)  !weight for high-symmetry q-points

      REAL*16,ALLOCATABLE :: enk(:), enkq(:,:)

      INTEGER :: ikx, iqx  ! index of k, q points
      REAL*8 :: distance

      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency
      INTEGER :: inu, nu
      REAL*8,ALLOCATABLE :: omega(:) ! real omega frequency
      REAL*8 :: omega_lower, omega_upper
      !COMPLEX*16,ALLOCATABLE: ommesh_phonon(:,:,:,:) ! phonon Matsubara frequency

      INTEGER :: i,j,l,m,n ! dummy index

      ! for MPI
      INTEGER rank,total,ierr,errorcode
      INTEGER pr_proc,numq_per,ikwp

      ! declare variables finished
      ! ###############################################################

      ! get input
      ! ###############################################################
      WRITE(*,*) 'do_EP_self starts:'

      INQUIRE(FILE = 'EP.in', EXIST = iffile)
      IF (iffile .EQV. .false.) THEN
          WRITE(*,*) 'input is needed for do_EP_self!'
          STOP
      ELSE
!           WRITE(*,*) 'Hi'
          OPEN(65, FILE = "EP.in")

          READ(65,*) T
          READ(65,*) mu
          READ(65,*) delta
          READ(65,*) nu
          READ(65,*) norb
          READ(65,*) eps_acustic
          READ(65,*) strength


          READ(65,*) numk
          READ(65,*) numomega, omega_lower, omega_upper
          READ(65,*) numq
          READ(65,*) numnv
          
!           WRITE(*,*) 'Ha'

!           !read in R-points
!           ALLOCATE(tran(nr,3))
!           READ(65,*) ((tran(i,j),j=1,3),i=1,nr)

!           !read in H(R)
!           ALLOCATE(ham(numr))
!           READ(65,*) (ham(i),i=1,numr)
          
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
          ALLOCATE(phonon_frequency(numq,numnv))
          READ(65,*) ((phonon_frequency(i,j),j=1,numnv),i=1,numq)

          WRITE(*,*) 'hello', phonon_frequency(1,1), phonon_frequency(1,4)
          
          !read in electron-phonon matrix g_v(k,q)
          ALLOCATE(Gv(numk,numq,numnv))
          READ(65,*) (((Gv(i,j,l),l=1,nu),j=1,numq),i=1,numk)
          
          !read in weight of each q-point w(q)
          ALLOCATE(weights(numq))
          READ(65,*) (weights(i),i=1,numq)

          CLOSE(65)
!           WRITE(*,*) 'Hello'
      ENDIF !IF iffile

      ! get input finished
      ! ###############################################################
      
      ! allocate variables
      ! ###############################################################
      ALLOCATE(sigi(numk))
      ALLOCATE(tot_sigi(numk))
      ! allocate variables finished


      ! calculate omega
      ! DO i = 1, numomega
      !  omega(i) = omega_lower + FLOAT(i-1)*(omega_upper-omega_lower)/FLOAT(numomega-1)

      ! ENDDO

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
      tot_sigi = 0.0

      WRITE(*,*) 'Rank: ', rank,  'pr_proc: ', pr_proc, ' numq_per: ', numq_per


      ! loop over k-points
      DO ikx = 1, numk

         ! loop over high-symmetry qpoints
         DO iqx = 1, numq

            ikwp = iqx-rank*pr_proc

            !IF(rank .EQ. 0) WRITE(*,*) 'strange: ', ikwp

            !parallel over different q points 
            IF (ikwp .GT. 0 .AND. ikwp .LE. numq_per) THEN 

               ! get epsilon_k. Fermi level is shifted.
               epsk = enk(ikx) - mu
               ! get epsilon_{k+q}. Fermi level is shifted.
               epskq = enkq(ikx,iqx) - mu

               DO inu = 1,nu
                
               !IF( rank .EQ. 0) WRITE(*,*) 'eps: ', wqnv, eps_acustic      

                  !gv(k,q)
                  temp1 = ABS(Gv(ikx,iqx,inu))

                  !phonon_frequency
                  wqnv = phonon_frequency(iqx,inu)
                  IF(rank .EQ. 0) WRITE(*,*) 'UUU: ', wqnv

                  IF( wqnv .GT. eps_acustic) THEN 

                    !Bose_Einstein occupancy
                    nb = 1.0/(EXP(wqnv/T)-1.0)

                    !Fermi-Dirac occupancy
                    nf = 1.0/((EXP(epskq/T)+1.0))
                                          
                    !(nb+nf)/(iomesh-ep+omega)                         
                    temp2 = (nb+nf)/(epsk+(0.0,1.0)*delta-epskq+wqnv)

                    !(nb+1-nf)/(iomesh-ep-omega)
                    temp3 = (nb+1.0-nf)/(epsk+(0.0,1.0)*delta-epskq-wqnv)

                    sigi(ikx)=sigi(ikx)+weights(iqx)*(temp1**2)*(temp2+temp3)

                            IF (n .EQ. 1 .AND. ABS(AIMAG(temp2)).LT.1.0) THEN
                             WRITE(*,*) 'hi: ',ikx,iqx,weights(iqx),inu
                             WRITE(*,*) 'why:', temp1, temp2, temp3, sigi(ikx)
                             WRITE(*,*) 'hii: ', omega(n), epskq, wqnv
                            ENDIF
  
                  ENDIF! wqnv

               ENDDO ! inu

            ENDIF ! ikwp

         ENDDO ! iqx

         ! use MPI to do the sum over q

         CALL MPI_REDUCE(sigi(ikx),tot_sigi(ikx),1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_Comm_World,ierr)

      ENDDO!ikx

      ! output
      IF (rank .EQ. 0) THEN

         OPEN(66, FILE = "EP.out")             
      
            WRITE(66, '(A,I4)') "number of k points: ", numk
            WRITE(66, '(A)') "K-points are: "
            DO i = 1, numk
               WRITE(66,'(3F24.16)') kpoints(i,1), kpoints(i,2), kpoints(i,3)
            ENDDO

         CLOSE(66)

         !write down spectral function

         distance = 0.0

         OPEN(67, FILE= "self-energy.out")
  
            DO ikx = 1, numk

               !WRITE(67,'(F5.2, F10.6, I, I, F24.12)') distance, omega(n), ikx, n, spectral(ikx,n)
               WRITE(67,'(2F10.6, 3F24.16)') distance, enk(ikx)-mu, REAL(tot_sigi(ikx)), AIMAG(tot_sigi(ikx))
               !WRITE(67,*) distance, omega(n), spectral(ikx,n), REAL(tot_sigi(ikx,n)), AIMAG(tot_sigi(ikx,n))

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
      DEALLOCATE(weights)
      DEALLOCATE(sigi)
      DEALLOCATE(tot_sigi)


      END PROGRAM do_EP_self
