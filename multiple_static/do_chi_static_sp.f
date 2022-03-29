      SUBROUTINE do_chi_static_sp(numq,q,T,Nc,mu,nk,nr,tran,ham,sigi_up,
     1sigi_down,chi_real,rank,total,norb,ncor_orb)
      !i\nu_m = 0, multiple orbital (norb*norb)

      !IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL*8, PARAMETER :: pi = 3.14159265358979
      REAL*8, PARAMETER :: smalleps = 1.e-5

      INTEGER :: numq !number of q points
      INTEGER :: iq !index of q
      REAL*8,DIMENSION(numq,3) :: q !q point
      REAL*8 :: qx, qy, qz !component of q point

      REAL*8 :: T !temperature, beta = 1/T
      INTEGER :: Nc !the cutoff for fermionic Matsubara frequency

      INTEGER :: nk(3), numk   ! number of k-point

      INTEGER :: c1, c2, c3
      INTEGER :: ikx, iky, ikz ! k index
      REAL*8  :: kx, ky, kz    ! k-point

      INTEGER :: i, i_shift    ! dummy index

      INTEGER :: n, n_shift    ! index of fermionic Matsubara frequency

      REAL*8,ALLOCATABLE :: ommesh(:)  ! fermionic Matsubara frequency

      INTEGER :: nr, r, norb, ncor_orb                 ! number of "r" points, dimension of ham
      INTEGER,DIMENSION(nr,3) :: tran  ! R
      COMPLEX*16,DIMENSION(nr,norb,norb) :: ham  ! H(R)

      COMPLEX*16,DIMENSION(norb,norb) :: hk, hkq ! noninteracting Hamiltonian at k and at k+q

      REAL*8,DIMENSION(numq) :: chi_real(numq), tot_chi_real(numq) ! irreducible susceptibility
      REAL*8,DIMENSION(numq) :: chi_real_up(numq), chi_real_down(numq)
      REAL*8,DIMENSION(numq) :: tot_chi_real_up(numq)
      REAL*8,DIMENSION(numq) :: tot_chi_real_down(numq)

      REAL*8 :: mu        ! chemical potential

      REAL*8, PARAMETER :: constant=0.5  ! for each spin channel
      COMPLEX*16,DIMENSION(ncor_orb,2*Nc) :: sigi_up, sigi_down ! self-energy
      ! REAL*8,DIMENSION(2*Nc) :: sigi_real,sigi_complex
      INTEGER,DIMENSION(2*Nc) :: w

      COMPLEX*16 :: trace_up,trace_down
      COMPLEX*16,ALLOCATABLE :: term_up(:,:),term_down(:,:)
      COMPLEX*16,DIMENSION(norb,norb) :: term1_up,term1_down
      ! COMPLEX*16,DIMENSION(norb,norb) :: part,term1_copy
      COMPLEX*16,DIMENSION(norb,norb) :: term2_up,term2_down
      COMPLEX*16,DIMENSION(norb,norb) :: identity,temp
      COMPLEX*16,DIMENSION(norb,norb) :: identity_sigma_up
      COMPLEX*16,DIMENSION(norb,norb) :: identity_sigma_down

      REAL*8,DIMENSION(numq) :: chi_imag, tot_chi_imag
      REAL*8,DIMENSION(numq) :: chi_imag_up, chi_imag_down
      REAL*8,DIMENSION(numq) :: tot_chi_imag_up, tot_chi_imag_down

      !for lapack and mpi
      INTEGER rank,total,ierr,errorcode
      INTEGER pr_proc,numk_per,ikp
      INTEGER :: info
      INTEGER,DIMENSION(norb) :: ipiv
      COMPLEX*16,DIMENSION(norb) :: work

      !real :: start_time, stop_time

      ALLOCATE(ommesh(2*Nc))

      !CALL cpu_time(start_time)

      !check symmetry
      ! DO i = -Nc,Nc-1
      ! 	i_shift = i + Nc + 1
      ! 	w(i_shift) = i
      ! 	sigi_real(i_shift) = REAL(sigi(i_shift))
      ! 	sigi_complex(i_shift) = AIMAG(sigi(i_shift))
      ! ENDDO

      ! set up fermionic Matsubara frequency
      DO i = -Nc, Nc-1
          i_shift = i + Nc + 1
          ommesh(i_shift) = (2*i+1)*pi*T
      ENDDO

      ! get the total number k-points
      numk = nk(1) * nk(2) * nk(3)

      DO i=0,11,1
          IF (i==rank) THEN
              WRITE(*,*),'Hello, World! I am process ',
     1        rank,' of ',total, '.'
          ENDIF
          CALL MPI_BARRIER(MPI_Comm_World,ierr)
      ENDDO

      ! pr_proc is the number of kpoints per CPU
      ! numk is equal to pr_proc except for the last CPU (highest rank)
      ! if nk/total is not an integer, numk = nk-rank*pr_proc
      pr_proc = FLOOR(numk/DBLE(total)+0.999)
      numk_per = pr_proc
      IF ((rank+1)*pr_proc .GT. numk) THEN
          IF (numk-rank*pr_proc .GT. 0) THEN
              numk_per = numk-rank*pr_proc
          ELSE
              numk_per = 0
          ENDIF
      ENDIF

!####################################################################################################

      ! initialize the quantities to be calculated
      chi_real_up = 0.0
      chi_real_down = 0.0
      chi_real = 0.0
      chi_imag_up = 0.0
      chi_imag_down = 0.0
      chi_imag = 0.0
      identity = 0.0
      identity_sigma_up = 0.0
      identity_sigma_down = 0.0

      DO i = 1, norb
          identity(i,i) = 1.0
      ENDDO

      DO iq = 1, numq
          qx = q(iq, 1)
          qy = q(iq, 2)
          qz = q(iq, 3)

      ! GG_real = 0.0
      ! GG_complex = 0.0

          DO ikx = 1, nk(1)
              DO iky = 1, nk(2)
                  DO ikz = 1, nk(3)

                      ikp = ikz+(iky-1)*nk(3)+(ikx-1)*nk(2)*nk(3)
     1                -rank*pr_proc

                      IF (ikp .GT. 0 .AND. ikp .LE. numk_per) THEN

                          kx = 2.0*pi*FLOAT(ikx)/FLOAT(nk(1))
                          ky = 2.0*pi*FLOAT(iky)/FLOAT(nk(2))
                          kz = 2.0*pi*FLOAT(ikz)/FLOAT(nk(3))

                          ! create H(k)...
                          hk = 0.0
                          DO r = 1,nr
                              temp=ham(r,:,:)
                              hk=hk+temp*EXP((0.0,1.0)*
     1                        (kx*tran(r,1)+ky*tran(r,2)+kz*tran(r,3)))
                          ENDDO !DO r
                          !WRITE(*,*) 'hk', hk

                          ! create H(k+q)
                          hkq = 0.0
                          DO r = 1,nr
                              temp=ham(r,:,:)
                              hkq=hkq+temp*EXP((0.0,1.0)*((kx+qx)*
     1                        tran(r,1)+(ky+qy)*tran(r,2)+
     1                        (kz+qz)*tran(r,3)))
                          ENDDO !DO r
                          !WRITE(*,*) 'hkq', hkq

                          term1_up = 0.0
                          term2_up = 0.0

                          term1_down = 0.0
                          term2_down = 0.0

                          DO n = 1, 2*Nc
                              !n_shift = n+(Nc+1)
                              n_shift = n

                              !term1 = 1.0/(((0.0,1.0)*ommesh(n_shift)+mu)*identity - hk - sigi(n_shift)) !G(k)
                              !term2 = 1.0/(((0.0,1.0)*ommesh(n_shift)+mu)*identity - hkq - sigi(n_shift)) !G(k+q)
                              DO i = 1, ncor_orb
                                  identity_sigma_up(i,i) = 
     1                            sigi_up(i,n_shift)
                                  identity_sigma_down(i,i) = 
     1                            sigi_down(i,n_shift)
                              ENDDO

                              term1_up = ((0.0,1.0)*ommesh(n_shift)
     1                        +mu)*identity - hk - identity_sigma_up !G(k)
                              term2_up = ((0.0,1.0)*ommesh(n_shift)
     1                        +mu)*identity - hkq - identity_sigma_up !G(k+q) vm=0

                              term1_down = ((0.0,1.0)*ommesh(n_shift)
     1                        +mu)*identity - hk - identity_sigma_down !G(k)
                              term2_down = ((0.0,1.0)*ommesh(n_shift)
     1                        +mu)*identity - hkq - identity_sigma_down !G(k+q) vm=0

                              !invert G(k)
                              !part = term1
                              !WRITE(*,*) 'term1 original', term1
                              CALL ZGETRF(norb, norb, term1_up, 
     1                        norb, ipiv, info)
                              !WRITE(*,*) 'term1 now', term1
                              !WRITE(*,*) 'info for F', info
                              CALL ZGETRI(norb, term1_up, norb, 
     1                        ipiv, work, norb, info)
                              !WRITE(*,*) 'should be one', matmul(term1,part)
                              !WRITE(*,*) 'info for I', info

                              CALL ZGETRF(norb, norb, term1_down, 
     1                        norb, ipiv, info)
                              CALL ZGETRI(norb, term1_down, norb, 
     1                        ipiv, work, norb, info)

                              !invert G(k+q)
                              !part = term2
                              CALL ZGETRF(norb, norb, term2_up, 
     1                        norb, ipiv, info)
                              !WRITE(*,*) 'info for F', info
                              CALL ZGETRI(norb, term2_up, norb, 
     1                        ipiv, work, norb, info)
                              !WRITE(*,*) 'should be one', matmul(term2,part)
                              !WRITE(*,*) 'info for I', info

                              CALL ZGETRF(norb, norb, term2_down, 
     1                        norb, ipiv, info)
                              CALL ZGETRI(norb, term2_down, norb, 
     1                        ipiv, work, norb, info)

                              CALL KPRODUCT(term1_up, term2_up,
     1                        term_up)
                              CALL KPRODUCT(term1_down, term2_down,
     1                        term_down)

                              trace_up = 0.0
                              trace_down = 0.0

                              DO i = 1, norb*norb
                                  trace_up = trace_up + term_up(i,i)
                                  trace_down = trace_down + 
     1                            term_down(i,i)
                              ENDDO

                              chi_real_up(iq) = chi_real_up(iq) 
     1                        + REAL(trace_up)
                              chi_imag_up(iq) = chi_imag_up(iq) 
     1                        + AIMAG(trace_up)

                              chi_real_down(iq) = chi_real_down(iq) 
     1                        + REAL(trace_down)
                              chi_imag_down(iq) = chi_imag_down(iq) 
     1                        + AIMAG(trace_down)

                              chi_real(iq) = chi_real_up(iq)*constant 
     1                        + chi_real_down(iq)*constant
                              chi_imag(iq) = chi_imag_up(iq)*constant 
     1                        + chi_imag_down(iq)*constant

                          ENDDO !DO n
                      ENDIF !DO ikp
                  ENDDO !DO ikz
              ENDDO !DO iky
          ENDDO !DO ikx
      ENDDO !iq

      !DO i = 1, 60
          !GG_real(i) = GG_real(i)/numk
          !GG_complex(i) = GG_complex(i)/numk
      !ENDDO

      !CALL cpu_time(stop_time)

      !WRITE(*,*) stop_time - start_time, "seconds"

!####################################################################################################
      CALL MPI_REDUCE(chi_real,tot_chi_real,numq,MPI_DOUBLE_PRECISION,
     1MPI_SUM,0,MPI_Comm_World,ierr)
      CALL MPI_REDUCE(chi_imag,tot_chi_imag,numq,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(chi_real_up,tot_chi_real_up,numq,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)
      CALL MPI_REDUCE(chi_real_down,tot_chi_real_down,numq,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      CALL MPI_REDUCE(chi_imag_up,tot_chi_imag_up,numq,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)
      CALL MPI_REDUCE(chi_imag_down,tot_chi_imag_down,numq,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      ! output
      IF (rank .EQ. 0) THEN
          DO i = 1, numq
              chi_real(i) = (-T/numk)*tot_chi_real(i)
              chi_imag(i) = (-T/numk)*tot_chi_imag(i)

              chi_real_up(i) = (-T/numk)*tot_chi_real_up(i)
              chi_imag_up(i) = (-T/numk)*tot_chi_imag_up(i)

              chi_real_down(i) = (-T/numk)*tot_chi_real_down(i)
              chi_imag_down(i) = (-T/numk)*tot_chi_imag_down(i)
          ENDDO

          DO i=1,numq
              WRITE(*,*) i, chi_real(i), chi_imag(i), chi_real_up(i),
     1        chi_imag_up(i), chi_real_down(i), chi_imag_down(i)
          ENDDO

      ENDIF !IF rank
      ! chi_real = (-T/numk)*tot_chi_real

      ! open(20, file='sigi.txt')

      ! close(20)

      DEALLOCATE(ommesh)
      RETURN

      END SUBROUTINE



