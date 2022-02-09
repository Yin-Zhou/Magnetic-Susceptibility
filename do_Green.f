      SUBROUTINE do_Green(nk,T,Nc,mu,ommesh,nr,tran,ham,sigi,Green,
     1           rank, total)
      !This function calculates the interacting Green function

      INCLUDE 'mpif.h'
      !IMPLICIT NONE

      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: smalleps = 1.e-5 

      INTEGER :: nk(3)  ! number of k-point
      REAL*8  :: T      ! temperature
      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency
      REAL*8  :: mu     ! chemical potential  
      REAL*8,DIMENSION(2*Nc) :: ommesh ! fermionic Matsubara frequency

      INTEGER :: nr                    ! number of "r" points
      INTEGER,DIMENSION(nr,3) :: tran  ! R
      COMPLEX*16,DIMENSION(nr) :: ham  ! H(R) 

      COMPLEX*16,DIMENSION(2*Nc) :: sigi ! self-energy
      COMPLEX*16,DIMENSION(Nc,nk(1),nk(2),nk(3)) :: Green !Green function

      INTEGER rank,total   ! for LAPACK and MPI

      ! local variables

      INTEGER :: ikx, iky, ikz   ! index of k
      INTEGER :: n     ! index of fermionic Matsubara frequency
      REAL*8  :: kx, ky, kz      ! component of k point
      INTEGER :: numk, numkw
      INTEGER :: r               ! number of "r" points
      COMPLEX*16 :: hk           ! noninteracting Hamiltonian at k and at k+q

      COMPLEX*16,ALLOCATABLE :: tot_Green(:,:,:,:) 

      ! for LAPACK and MPI
      INTEGER ierr,errorcode
      INTEGER pr_proc,numwk,ikwp

      ALLOCATE(tot_Green(Nc, nk(1), nk(2), nk(3)))

      WRITE(*,*) 'entering do_Green'

      numk  = nk(1)*nk(2)*nk(3)
      numkw = nk(1)*nk(2)*nk(3)*Nc

      ! pr_proc is the number of kpoints per CPU
      ! numk is equal to pr_proc except for the last CPU (highest rank)
      ! if nk/total is not an integer, numk = nk-rank*pr_proc

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
     1   numkw_per

      ! initialize the quantities to be calculated
      Green = 0.0

      WRITE(*,*), 'hereNew 0'

      DO n = 1, Nc

         DO ikx = 1, nk(1) 

            DO iky = 1, nk(2)

               DO ikz = 1, nk(3)

                   ikwp = ikz+(iky-1)*nk(3)+(ikx-1)*nk(2)*nk(3)+
     1             (n-1)*nk(1)*nk(2)*nk(3)-rank*pr_proc

                   IF (ikwp .GT. 0 .AND. ikwp .LE. numkw_per) THEN
                     
                       kx = 2.0*pi*FLOAT(ikx)/FLOAT(nk(1))
                       ky = 2.0*pi*FLOAT(iky)/FLOAT(nk(2))
                       kz = 2.0*pi*FLOAT(ikz)/FLOAT(nk(3))

                       ! create H(k)...
                       hk = 0.0
                       DO r = 1,nr
                       !WRITE(*,*) 'tran: ', tran(r,1)
                           hk=hk+ham(r)*EXP(
     1               (0.0,1.0)*(kx*tran(r,1)+ky*tran(r,2)+kz*tran(r,3)))
                       ENDDO !DO r
                      
                       ! calculate G(k)
                       Green(n,ikx,iky,ikz) = 
     1               1.0/((0.0,1.0)*ommesh(n)+mu-hk-sigi(n))

                   ENDIF !DO ikwp
               ENDDO !DO ikz
            ENDDO !DO iky
         ENDDO !DO ikx

      ENDDO !n

      CALL MPI_ALLREDUCE(Green,tot_Green,Nc*nk(1)*nk(2)*nk(3),
     1 MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_Comm_World,ierr)
 
      Green = (-T/numk)*tot_Green 

      DEALLOCATE(tot_Green)                     

      RETURN 

      END SUBROUTINE