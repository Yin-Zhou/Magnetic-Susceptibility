      PROGRAM do_EPsuper
      INCLUDE 'mpif.h'
      !IMPLICIT NONE      
      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: smalleps = 1.e-5

      ! declare variables
      ! ###############################################################
      INTEGER :: nq(3)  ! number of q points
      INTEGER :: iqx,iqy,iqz     ! index of q points

      INTEGER :: nk(3) ! number of k-point

      COMPLEX*16,ALLOCATABLE :: Gv(:,:,:,:,:,:,:)
      COMPLEX*16,ALLOCATABLE :: sigi(:,:,:,:,:,:)
      COMPLEX*16,ALLOCATABLE :: tot_sigi(:,:,:,:,:,:)
      COMPLEX*16 :: temp1,temp2
      COMPLEX*16,DIMENSION(2,2) :: temp3
      COMPLEX*16,ALLOCATABLE :: GREEN(:,:,:,:,:,:)
      COMPLEX*16,ALLOCATABLE :: tot_GREEN(:,:,:,:,:,:)
      COMPLEX*16,ALLOCATABLE :: GREENAC(:,:,:,:,:,:)
      COMPLEX*16,ALLOCATABLE :: tot_GREENAC(:,:,:,:,:,:)
      REAL*8,ALLOCATABLE :: spectral(:,:,:,:)
      REAL*8,ALLOCATABLE :: tot_spectral(:,:,:,:)
      REAL*8 :: nb,nf,wqmiu,epskq
      REAL*8,ALLOCATABLE :: weights(:,:,:)
      REAL*8 :: weight
      REAL*8,ALLOCATABLE :: delta_k(:,:,:)

      COMPLEX*16,DIMENSION(2,2) :: tau0,tau1,tau2,tau3

      REAL*8  :: k(3)  ! k and k'
      INTEGER :: ikx, iky, ikz  ! index of k, k' points

      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency
      INTEGER :: Ncq, miu, imiu
      REAL*8,ALLOCATABLE :: ommesh(:)  ! fermionic Matsubara frequency
      REAL*8,ALLOCATABLE :: ommeshq(:)
      REAL*8,ALLOCATABLE :: OMEGA(:,:,:,:)
      REAL*8,ALLOCATABLE :: eps(:,:,:)

      INTEGER :: i,j,l,m,n,x,y ! dummy index
      INTEGER :: numk,numkw

      ! for MPI
      INTEGER rank,total,ierr,errorcode
      INTEGER pr_proc,numkw_per,ikwp

      ! declare variables finished
      ! ###############################################################

      ! get input
      ! ###############################################################
      WRITE(*,*) 'do_EPsuper starts:'

      INQUIRE(FILE = 'do_EPsuper.in', EXIST = iffile)
      IF (iffile .EQV. .false.) THEN
          WRITE(*,*) 'input is needed for do_EPsuper!'
          STOP
      ELSE
          OPEN(65, FILE = "do_EPsuper.in")

          READ(65,*) T
          READ(65,*) mu

          READ(65,*) (nk(i),i=1,3)
          READ(65,*) Nc

          READ(65,*) (nq(i),i=1,3)
          READ(65,*) Ncq

          READ(65,*) miu

          ALLOCATE(Gv(nk(1),nk(2),nk(3),nq(1),nq(2),nq(3),miu))
          READ(65,*) (((((((Gv(i,j,l,m,n,x,y),y=1,miu),x=1,nq(3)),
     1    n=1,nq(2)),m=1,nq(1)),l=1,nk(3)),j=1,nk(2)),i=1,nk(1))
          
          ALLOCATE(weights(nq(1),nq(2),nq(3)))
          READ(65,*) (((weights(i,j,x),x=1,nq(3)),j=1,
     1    nq(2)),i=1,nq(1))

          ALLOCATE(eps(nk(1)+nq(1),nk(2)+nq(2),nk(3)+nq(3)))
          READ(65,*) (((eps(i,j,x),x=1,nk(3)+nq(3)),j=1,nk(2)+nq(2)),
     1    i=1,nk(1)+nq(1))

          ALLOCATE(OMEGA(nq(1),nq(2),nq(3),miu))
          READ(65,*) ((((OMEGA(i,j,x,y),y=1,miu),x=1,nq(3)),j=1,nq(2)),
     1    i=1,nq(1))

          ALLOCATE(delta_k(nk(1),nk(2),nk(3)))
          READ(65,*) (((delta_k(i,j,x),x=1,nk(3)),j=1,nk(2)),i=1,nk(1))         

          CLOSE(65)
      ENDIF !IF iffile

      ! get input finished
      ! ###############################################################
      
      ! allocate variables
      ! ###############################################################
      ALLOCATE(ommesh(2*Nc))
      ALLOCATE(ommeshq(2*Ncq))
      ALLOCATE(sigi(nk(1),nk(2),nk(3),2*Nc,2,2))
      ALLOCATE(tot_sigi(nk(1),nk(2),nk(3),2*Nc,2,2))
      ALLOCATE(GREEN(nk(1),nk(2),nk(3),2*Nc,2,2))
      ALLOCATE(tot_GREEN(nk(1),nk(2),nk(3),2*Nc,2,2))
      ALLOCATE(GREENAC(nk(1),nk(2),nk(3),2*Nc,2,2))
      ALLOCATE(tot_GREENAC(nk(1),nk(2),nk(3),2*Nc,2,2))
      ALLOCATE(spectral(nk(1),nk(2),nk(3),2*Nc))
      ALLOCATE(tot_spectral(nk(1),nk(2),nk(3),2*Nc))
      ! allocate variables finished
      ! ###############################################################

      tau0 = reshape((/1.0,0.0,0.0,1.0/),shape(tau0))
      tau1 = reshape((/0.0,1.0,1.0,0.0/),shape(tau1))
      tau2 = reshape((/0.0,(0.0,1.0),(0.0,-1.0),0.0/),shape(tau2))
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

      numk  = nk(1)*nk(2)*nk(3)
      numkw = nq(1)*nq(2)*nq(3)*(2*Ncq)
 
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

      sigi = 0.0
      GREEN = 0.0
      GREENAC = 0.0

      WRITE(*,*) 'Rank: ', rank,  'pr_proc: ', pr_proc, ' numkw_per: ',
     1numkw_per

      DO n = -Nc, Nc-1
          n_shift = n + (Nc+1)
          DO ikx = 1, nk(1)
              DO iky = 1, nk(2)
                  DO ikz = 1, nk(3)

                      ikwp = ikz+(iky-1)*nk(3)+(ikx-1)*nk(2)*nk(3)+
     1                (n_shift-1)*nk(1)*nk(2)*nk(3)-rank*pr_proc

                      ! parallel over different k points and omega
                      IF (ikwp .GT. 0 .AND. ikwp .LE. numkw_per) THEN     

                          DO m = -Ncq, Ncq-1

                              m_shift = m + (Ncq+1)

                              DO iqx = 1, nq(1)
                                  DO iqy = 1, nq(2)
                                      DO iqz = 1, nq(3)

                                          weight=weights(iqx,iqy,iqz)

                                          DO imiu = 1,miu
                                          
                                          ! gv(k,q)
                                          temp1 = Gv(ikx,iky,ikz,iqx,
     1                                    iqy,iqz,imiu)

                                          !w q miu
                                          wqmiu = OMEGA(iqx,iqy,
     1                                    iqz,imiu)

                                          !eps k+q
                                          epsk = eps(ikx,iky,ikz)
                                          
                                          ! 2omega/(iNcq**2-omega**2)
                                          temp2 = (2*wqmiu)/(((0.0,1.0)
     1                                    *ommeshq(m_shift)**2)-
     1                                    wqmiu**2)

                                          ! iomesh tau0 + ep tau3 - delta_k tau1
                                          temp3 = ((0.0,1.0)*
     1                                    ommesh(n_shift)*tau0+
     1                                    epsk*tau3-
     1                                    delta_k(ikx,iky,ikz)*tau1)/
     1                                    ((0.0,1.0)*
     1                                    ommesh(n_shift)**2+
     1                                    SQRT(epsk**2+
     1                                    delta_k(ikx,iky,ikz)**2))

                                          sigi(ikx,iky,ikz,n_shift,:,:)=
     1                                    sigi(ikx,iky,ikz,n_shift,:,:)+
     1                                    weight*(temp1**2)*temp2*
     1                                    temp3

                                          ENDDO ! imiu
                                      ENDDO ! iqz
                                  ENDDO ! iqy
                              ENDDO ! iqx
                          ENDDO ! m
                      ENDIF !ikwp

                      sigi(ikx,iky,ikz,n_shift,:,:) = -T*
     1                sigi(ikx,iky,ikz,n_shift,:,:)

                      GREEN(ikx,iky,ikz,n_shift,:,:) = 1/((0.0,1.0)*
     1                ommesh(n_shift)*tau0-eps(ikx,iky,ikz)*tau3-
     1                sigi(ikx,iky,ikz,n_shift,:,:))

                      GREENAC(ikx,iky,ikz,n_shift,:,:) = 1/(
     1                ommesh(n_shift)*tau0+(0.0,1.0)*smalleps*tau0-
     1                eps(ikx,iky,ikz)*tau3-
     1                sigi(ikx,iky,ikz,n_shift,:,:))

                      spectral(ikx,iky,ikz,n_shift) = 
     1                -AIMAG(GREENAC(ikx,iky,ikz,n_shift,1,1))/pi        

                  ENDDO ! ikz
              ENDDO ! iky
          ENDDO ! ikx
      ENDDO ! n

      CALL MPI_REDUCE(sigi,tot_sigi,nk(1)*nk(2)*nk(3)*2*Nc*4,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(GREEN,tot_GREEN,nk(1)*nk(2)*nk(3)*2*Nc*4,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(GREENAC,tot_GREENAC,nk(1)*nk(2)*nk(3)*2*Nc*4,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(spectral,tot_spectral,nk(1)*nk(2)*nk(3)*2*Nc,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr) 

      ! output
      IF (rank .EQ. 0) THEN

      OPEN(65, FILE = "DMFT_Eliashberg.out")             
      
      WRITE(65,'(4I3,1F24.16)') 2*Nc,nk(1),nk(2),nk(3),
     1 REAL(tot_spectral(nk(1),nk(2),nk(3),2*Nc))
      CLOSE(65)

      ENDIF ! rank     

      CALL MPI_Finalize(ierr)      

      DEALLOCATE(ommesh)
      DEALLOCATE(ommeshq)
      DEALLOCATE(Gv)
      DEALLOCATE(weights)
      DEALLOCATE(sigi)
      DEALLOCATE(tot_sigi)
      DEALLOCATE(eps)
      DEALLOCATE(OMEGA)
      DEALLOCATE(delta_k)
      DEALLOCATE(GREEN)
      DEALLOCATE(tot_GREEN)
      DEALLOCATE(GREENAC)
      DEALLOCATE(tot_GREENAC)
      DEALLOCATE(spectral)
      DEALLOCATE(tot_spectral)

      END PROGRAM do_EPsuper