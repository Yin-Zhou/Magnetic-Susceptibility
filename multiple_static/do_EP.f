      PROGRAM do_EP
      INCLUDE 'mpif.h'
      !IMPLICIT NONE      
      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: smalleps = 1.e-5

      ! declare variables
      ! ###############################################################
      INTEGER :: nq(3)  ! number of q points
      INTEGER :: iqx,iqy,iqz     ! index of q points

      INTEGER :: nk,Nomega ! number of k-point

      COMPLEX*16,ALLOCATABLE :: Gv(:,:,:,:,:)
      COMPLEX*16,ALLOCATABLE :: sigi(:,:)
      COMPLEX*16,ALLOCATABLE :: tot_sigi(:,:)
      COMPLEX*16 :: temp1,temp2,temp3
      COMPLEX*16,ALLOCATABLE :: GREEN(:,:)
      COMPLEX*16,ALLOCATABLE :: tot_GREEN(:,:)
      COMPLEX*16,ALLOCATABLE :: GREENAC(:,:)
      COMPLEX*16,ALLOCATABLE :: tot_GREENAC(:,:)
      REAL*8,ALLOCATABLE :: spectral(:,:)
      REAL*8,ALLOCATABLE :: tot_spectral(:,:)
      REAL*8 :: nb,nf,wqmiu,epskq
      REAL*8,ALLOCATABLE :: weights(:,:,:)
      REAL*8 :: weight

C       REAL*8  :: nk  ! k and k'
      INTEGER :: ikx, iky, ikz  ! index of k
      REAL*8,ALLOCATABLE :: kpoints(:)

      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency
C       INTEGER :: Ncq
      INTEGER :: nu, inu
      REAL*8,ALLOCATABLE :: ommesh(:)  ! fermionic Matsubara frequency
      REAL*8,ALLOCATABLE :: OMEGA(:,:,:,:)
      REAL*8,ALLOCATABLE :: eps(:)

      INTEGER :: i,j,l,m,n,x,y ! dummy index
      INTEGER :: numk,numkw

      ! for MPI
      INTEGER rank,total,ierr,errorcode
      INTEGER pr_proc,numkw_per,ikwp

      ! declare variables finished
      ! ###############################################################

      ! get input
      ! ###############################################################
      WRITE(*,*) 'do_EP starts:'

      INQUIRE(FILE = 'do_EP.in', EXIST = iffile)
      IF (iffile .EQV. .false.) THEN
          WRITE(*,*) 'input is needed for do_EP!'
          STOP
      ELSE
          OPEN(65, FILE = "do_EP.in")

          READ(65,*) T
          READ(65,*) mu

          READ(65,*) nk
          READ(65,*) Nomega

          READ(65,*) (nq(i),i=1,3)
C           READ(65,*) Ncq

          READ(65,*) nu !number of v

          ALLOCATE(kpoints(2*nk))
          READ(65,*) (kpoints(i),i=1,2*nk)

          ALLOCATE(ommesh(Nomega))
          READ(65,*) (ommesh(i),i=1,Nomega)

          ALLOCATE(Gv(nk,nq(1),nq(2),nq(3),nu))
          READ(65,*) (((((Gv(i,j,l,m,n),n=1,nu),m=1,nq(3)),l=1,nq(2)),
     1    j=1,nq(1)),i=1,nk)
          
          ALLOCATE(weights(nq(1),nq(2),nq(3)))
          READ(65,*) (((weights(i,j,x),x=1,nq(3)),j=1,
     1    nq(2)),i=1,nq(1))

          ALLOCATE(eps(nk+nq(1)+nq(2)+nq(3)))
          READ(65,*) (eps(i),x=1,nk+nq(1)+nq(2)+nq(3))

          ALLOCATE(OMEGA(nq(1),nq(2),nq(3),nu))
          READ(65,*) ((((OMEGA(i,j,x,y),y=1,nu),x=1,nq(3)),j=1,nq(2)),
     1    i=1,nq(1))         

          CLOSE(65)
      ENDIF !IF iffile

      ! get input finished
      ! ###############################################################
      
      ! allocate variables
      ! ###############################################################
      ALLOCATE(sigi(nk,Nomega))
      ALLOCATE(tot_sigi(nk,Nomega))
      ALLOCATE(GREEN(nk,Nomega))
      ALLOCATE(tot_GREEN(nk,Nomega))
      ALLOCATE(GREENAC(nk,Nomega))
      ALLOCATE(tot_GREENAC(nk,Nomega))
      ALLOCATE(spectral(nk,Nomega))
      ALLOCATE(tot_spectral(nk,Nomega))
      ! allocate variables finished
      ! ###############################################################

C       DO i = -Nc, Nc-1
C           i_shift = i + (Nc+1)
C           ommesh(i_shift) = (2*i+1)*pi*T
C       ENDDO

      ! launch MPI
      CALL MPI_init(ierr)

      IF (ierr .NE. MPI_SUCCESS) THEN
          WRITE(*,*), 'Error starting MPI program. Terminating.'
          CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      ENDIF !IF ierr

      CALL MPI_Comm_size(MPI_Comm_World, total, ierr)
      CALL MPI_Comm_rank(MPI_Comm_World, rank, ierr)

      numk  = nk
      numkw = nq(1)*nq(2)*nq(3)

      ! pr_proc is the number of kpoints per CPU
      ! numk_per is equal to pr_proc except for the last CPU (highest
      ! rank)
      ! if nk/total is not an integer, numk_per = numk-rank*pr_proc

      pr_proc = FLOOR(numk/DBLE(total)+0.999)

      numkw_per = pr_proc
      IF ((rank+1)*pr_proc .GT. numk) THEN
          IF (numk-rank*pr_proc .GT. 0) THEN
              numk_per = numk-rank*pr_proc
          ELSE
              numk_per = 0
          ENDIF
      ENDIF

      sigi = 0.0
      GREEN = 0.0
      GREENAC = 0.0

      WRITE(*,*) 'Rank: ', rank,  'pr_proc: ', pr_proc, ' numkw_per: ',
     1numkw_per

      DO n = 1, Nomega
          n_shift = n
          DO ikx = 1, nk
              ikwp = (n_shift-1)*nk-rank*pr_proc

              ! parallel over different k points and omega
              IF (ikwp .GT. 0 .AND. ikwp .LE. numkw_per) THEN     

C                           DO m = -Ncq, Ncq-1

C                               m_shift = m + (Ncq+1)

                  DO iqx = 1, nq(1)
                      DO iqy = 1, nq(2)
                          DO iqz = 1, nq(3)

                              weight=weights(iqx,iqy,iqz)

                              DO inu = 1,nu
                                          
                                  ! gv(k,q)
                                  temp1 = Gv(ikx,iqx,iqy,iqz,inu)

                                  !w q miu
                                  wqmiu = OMEGA(iqx,iqy,iqz,inu)

                                  !eps k+q
                                  epskq = eps(ikx+iqx+iqy+iqz)

                                  !nb
                                  nb = 1/(EXP(wqmiu/T)-1)

                                  !nf
                                  nf = 1/((EXP((epskq-mu)/T)+1))
                                          
                                  ! (nb+nf)/(iomesh-ep+omega)
                                  temp2 = (nb+nf)/((0.0,1.0)*
     1                            ommesh(n_shift)-epskq+wqmiu)

                                  ! (nb+1-nf)/(iomesh-ep-omega)
                                  temp3 = (nb+1-nf)/((0.0,1.0)*
     1                            ommesh(n_shift)-epskq-wqmiu)

                                  sigi(ikx,n_shift)=
     1                            sigi(ikx,n_shift)+
     1                            weight*(temp1**2)*(temp2+temp3)

                              ENDDO ! inu
                          ENDDO ! iqz
                      ENDDO ! iqy
                  ENDDO ! iqx
              ENDIF !ikwp

C                       GREEN(ikx,iky,ikz,n_shift) = 1/((0.0,1.0)*
C      1                ommesh(n_shift)-eps(ikx,iky,ikz)-
C      1                sigi(ikx,iky,ikz,n_shift))

              GREENAC(ikx,n_shift) = 1/(ommesh(n_shift)+
     1        (0.0,1.0)*smalleps-eps(ikx)-
     1        sigi(ikx,n_shift))

              spectral(ikx,n_shift) = 
     1        -AIMAG(GREENAC(ikx,n_shift))/pi        

          ENDDO ! ikx
      ENDDO ! n

      CALL MPI_REDUCE(sigi,tot_sigi,nk*Nomega,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(GREEN,tot_GREEN,nk*Nomega,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(GREENAC,tot_GREENAC,nk*Nomega,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(spectral,tot_spectral,nk*Nomega,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr) 

      ! output
      IF (rank .EQ. 0) THEN

      OPEN(65, FILE = "DMFT.out")             
      
      WRITE(65,'(2I3,1F24.16)') nk,Nomega,
     1 REAL(tot_spectral(nk,Nomega))
      CLOSE(65)

      ENDIF ! rank     

      CALL MPI_Finalize(ierr)      

      DEALLOCATE(kpoints)
      DEALLOCATE(ommesh)
      DEALLOCATE(Gv)
      DEALLOCATE(weights)
      DEALLOCATE(sigi)
      DEALLOCATE(tot_sigi)
      DEALLOCATE(eps)
      DEALLOCATE(OMEGA)
      DEALLOCATE(GREEN)
      DEALLOCATE(tot_GREEN)
      DEALLOCATE(GREENAC)
      DEALLOCATE(tot_GREENAC)
      DEALLOCATE(spectral)
      DEALLOCATE(tot_spectral)

      END PROGRAM do_EP

