      SUBROUTINE do_chi(nq,T,Nc,mu,nk,ommesh,nr,tran,ham,sigi,
     1                  chi,rank,total)
      !This function calculates irreducible susceptibility chi
      !at a given q and inv_m
      !The information of "q" is hidden in "hkq".

      INCLUDE 'mpif.h'
      !IMPLICIT NONE

      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: smalleps = 1.E-5

      INTEGER :: nq(3)         ! number of q points
      INTEGER :: iqx,iqy,iqz   ! index of q
      REAL*8  :: qx, qy, qz    ! component of q point

      INTEGER :: m, m_shift    ! index of bosonic Matsubara frequency
                               ! here 'm' should be non-negative
      INTEGER :: numnu

      REAL*8  :: T      ! temperature
      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency

      INTEGER :: nk(3)  ! number of k-point
      INTEGER :: numk, numqw

      INTEGER :: ikx, iky, ikz ! index
      REAL*8  :: kx, ky, kz    ! k-point

      INTEGER :: i, i_shift    ! dummy index

      INTEGER :: n, n_shift    ! index of fermionic Matsubara frequency

      REAL*8,DIMENSION(2*Nc) :: ommesh ! fermionic Matsubara frequency

      INTEGER :: nr, r                 ! number of "r" points
      INTEGER,DIMENSION(nr,3) :: tran  ! R
      COMPLEX*16,DIMENSION(nr) :: ham  ! H(R)

      COMPLEX*16 :: hk, hkq ! noninteracting Hamiltonian at k and at k+q

      COMPLEX*16,DIMENSION(2*Nc,nq(1),nq(2),nq(3)) :: chi ! irreducible susceptibility 

      COMPLEX*16,ALLOCATABLE :: tot_chi(:,:,:,:)

      REAL*8 :: mu        ! chemical potential  
      REAL*8 :: fermi     ! Fermi-Dirac occupancy
      REAL*8 :: fermi_dev ! Derivative of Fermi-Dirac occupancy

      COMPLEX*16,DIMENSION(2*Nc) :: sigi ! self-energy

      COMPLEX*16 :: term1, term2

      ! for LAPACK and MPI
      INTEGER rank,total,ierr,errorcode
      INTEGER pr_proc,numqw_per,iqwp

      ALLOCATE(tot_chi(2*Nc, nk(1), nk(2), nk(3)))

      WRITE(*,*) 'entering do_chi'

      ! get the total number k-points
      numk = nk(1) * nk(2) * nk(3)
      
      numqw = nq(1)*nq(2)*nq(3)*(2*Nc)
      
      !WRITE(*,*),'Hello, World! I am process ',rank,' of ',total, '.'

      ! pr_proc is the number of kpoints per CPU
      ! numk_per is equal to pr_proc except for the last CPU (highest rank)
      ! if nk/total is not an integer, numk_per = numk-rank*pr_proc
      pr_proc = FLOOR(numqw/DBLE(total)+0.999)

      numqw_per = pr_proc
      IF ((rank+1)*pr_proc .GT. numqw) THEN
          IF (numqw-rank*pr_proc .GT. 0) THEN
              numqw_per = numqw-rank*pr_proc
          ELSE
              numqw_per = 0
          ENDIF
      ENDIF

      !debug
      WRITE(*,*), 'Rank: ',rank,'pr_proc: ',pr_proc,'numqw_per: ', 
     1numqw_per

      ! initialize the quantities to be calculated
      chi = 0.0

      WRITE(*,*), 'hereNew 0'

      numnu = 2*Nc -1

      DO m = 0, numnu 

          m_shift = m + 1

          DO iqx = 1, nq(1)
              DO iqy = 1, nq(2)
                  DO iqz = 1, nq(3)

                      qx = 2.0*pi*FLOAT(iqx)/FLOAT(nq(1))
                      qy = 2.0*pi*FLOAT(iqy)/FLOAT(nq(2))
                      qz = 2.0*pi*FLOAT(iqz)/FLOAT(nq(3))

                      iqwp = iqz+(iqy-1)*nq(3)+(iqx-1)*nq(2)*nq(3)+
     1                (m_shift-1)*nq(1)*nq(2)*nq(3)-rank*pr_proc

                      IF (iqwp .GT. 0 .AND. iqwp .LE. numqw_per) THEN

C                           IF (rank .eq. 0) WRITE(*,*) 'iqwp: ', iqwp

                          DO ikx = 1, nk(1)
                              DO iky = 1, nk(2)
                                  DO ikz = 1, nk(3)

                                      kx = 2.0*pi*
     1                                FLOAT(ikx)/FLOAT(nk(1))
                                      ky = 2.0*pi*
     1                                FLOAT(iky)/FLOAT(nk(2))
                                      kz = 2.0*pi*
     1                                FLOAT(ikz)/FLOAT(nk(3))


                                      ! create H(k)...
                                      hk = 0.0
                                      DO r = 1,nr
                                          hk=hk+ham(r)*EXP(
     1              (0.0,1.0)*(kx*tran(r,1)+ky*tran(r,2)+kz*tran(r,3)))
                                      ENDDO !DO r

                                      ! create H(k+q)
                                      hkq = 0.0
                                      DO r = 1,nr
                                          hkq = hkq+ham(r)*EXP(
     1(0.0,1.0)*((kx+qx)*tran(r,1)+(ky+qy)*tran(r,2)+(kz+qz)*tran(r,3)))
                                      ENDDO !DO r

                                      term1 = 0.0
                                      term2 = 0.0

                                      DO n = -Nc, Nc-1
                                          n_shift = n+(Nc+1)

                                          term1 = 1.0/((0.0,1.0)*
     1                         ommesh(n_shift+m)+mu-hkq-sigi(n_shift+m))

                                          term2 = 1.0/((0.0,1.0)*
     1                         ommesh(n_shift)+mu-hk-sigi(n_shift))

                               chi(m_shift,iqx,iqy,iqz) = 
     1                         chi(m_shift,iqx,iqy,iqz)+term1*term2
                                      ENDDO ! n
                                  ENDDO ! ikz
                              ENDDO ! iky
                          ENDDO ! ikx
                      ENDIF ! iqwp
                  ENDDO ! iqz 
              ENDDO ! iqy
          ENDDO ! iqx
      ENDDO !m

      CALL MPI_ALLREDUCE(chi,tot_chi,(numnu+1)*nq(1)*nq(2)*nq(3),
     1 MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_Comm_World,ierr)

      chi = (-T/numk)*tot_chi

      DEALLOCATE(tot_chi)

      RETURN 

      END SUBROUTINE