      SUBROUTINE do_chi_static(numq,q,T,Nc,mu,nk,nr,tran,ham,sigi,U_S,
     1 CHI_S,CHI_C,chi_real,rank,total,norb,ncor_orb)
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

      INTEGER :: i, i_shift   ! dummy index
      INTEGER :: l1, l2, l3, l4
      INTEGER :: d1, d2
      INTEGER :: c4
      INTEGER :: r_in, c_in, index_r, index_c

      INTEGER :: n, n_shift    ! index of fermionic Matsubara frequency

      REAL*8,ALLOCATABLE :: ommesh(:)  ! fermionic Matsubara frequency            

      INTEGER :: nr, r, norb, ncor_orb                 ! number of "r" points, dimension of ham
      INTEGER,DIMENSION(nr,3) :: tran  ! R
      COMPLEX*16,DIMENSION(nr,norb,norb) :: ham  ! H(R)

      COMPLEX*16,DIMENSION(norb,norb) :: hk, hkq ! noninteracting Hamiltonian at k and at k+q
      
      ! irreducible susceptibility
      ! chi tensor with 4 indices
      REAL*8,DIMENSION(norb,norb,norb,norb,numq) :: chi_real 
      REAL*8,DIMENSION(norb,norb,norb,norb,numq) :: tot_chi_real

      REAL*8,DIMENSION(norb,norb,norb,norb,numq) :: chi_imag 
      REAL*8,DIMENSION(norb,norb,norb,norb,numq) :: tot_chi_imag

      !create the matrix version of chi
      COMPLEX*16,DIMENSION(norb*norb,norb*norb,numq) :: CHI_0
      COMPLEX*16,DIMENSION(norb*norb,norb*norb) :: U_S
      COMPLEX*16,DIMENSION(norb*norb,norb*norb,numq) :: CHI_S
      COMPLEX*16,DIMENSION(norb*norb,norb*norb,numq) :: CHI_C
      !COMPLEX*16,DIMENSION(norb*norb,norb*norb,numq) :: tot_CHI_S
      COMPLEX*16,DIMENSION(norb*norb,norb*norb) :: part

      REAL*8 :: mu        ! chemical potential

      COMPLEX*16,DIMENSION(ncor_orb,2*Nc) :: sigi! self-energy
      ! REAL*8,DIMENSION(2*Nc) :: sigi_real,sigi_complex
      !INTEGER,DIMENSION(2*Nc) :: w

      !COMPLEX*16 :: trace
      COMPLEX*16,DIMENSION(norb,norb,norb,norb):: term
      ! COMPLEX*16,DIMENSION(norb,norb) :: part,term1_copy
      COMPLEX*16,DIMENSION(norb,norb) :: term1,term2,identity,temp
      COMPLEX*16,DIMENSION(norb,norb) :: identity_sigma 
      COMPLEX*16,DIMENSION(norb*norb,norb*norb) :: identity_S

      !for lapack and mpi
      INTEGER rank,total,ierr
      !INTEGER errorcode
      INTEGER pr_proc,numk_per,ikp
      INTEGER :: info
      INTEGER,DIMENSION(norb) :: ipiv
      COMPLEX*16,DIMENSION(norb) :: work

      INTEGER,DIMENSION(norb*norb) :: ipivS
      COMPLEX*16,DIMENSION(norb*norb) :: workS

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

      DO i=0,total,1
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
      
      !WRITE(*,*) 'Hello'
!####################################################################################################

      ! initialize the quantities to be calculated
      chi_real = 0.0
      chi_imag = 0.0
      identity = 0.0
      identity_sigma = 0.0

      CHI_0 = 0.0
      CHI_S = 0.0
      CHI_C = 0.0

      DO i = 1, norb
          identity(i,i) = 1.0
      ENDDO

      identity_S = 0.0
      DO i = 1, norb*norb
          identity_S(i,i) = 1.0
      ENDDO !do i

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
                          !WRITE(*,*) 'Hey'

                          term1 = 0.0
                          term2 = 0.0
                          !part = 0.0

                          DO n = 1, 2*Nc

                              !n_shift = n+(Nc+1)
                              n_shift = n

                              !term1 = 1.0/(((0.0,1.0)*ommesh(n_shift)+mu)*identity - hk - sigi(n_shift)) !G(k)
                              !term2 = 1.0/(((0.0,1.0)*ommesh(n_shift)+mu)*identity - hkq - sigi(n_shift)) !G(k+q)
                              DO i = 1, ncor_orb
                                  identity_sigma(i,i) = 
     1                            sigi(i,n_shift)
                              ENDDO

                              term1 = ((0.0,1.0)*ommesh(n_shift)
     1                        +mu)*identity - hk - identity_sigma !G(k)
                              term2 = ((0.0,1.0)*ommesh(n_shift)
     1                        +mu)*identity - hkq - identity_sigma !G(k+q) vm=0

                              !invert G(k)
                              !part = term1
                              !WRITE(*,*) 'term1 original', term1
                              CALL ZGETRF(norb, norb, term1, 
     1                        norb, ipiv, info)
                              !WRITE(*,*) 'term1 now', term1
                              !WRITE(*,*) 'info for F', info
                              CALL ZGETRI(norb, term1, norb, 
     1                        ipiv, work, norb, info)
                              !WRITE(*,*) 'should be one', matmul(term1,part)
                              !WRITE(*,*) 'info for I', info

                              !invert G(k+q)
                              !part = term1
                              !WRITE(*,*) 'term1 original', term1
                              CALL ZGETRF(norb, norb, term2, 
     1                        norb, ipiv, info)
                              !WRITE(*,*) 'term1 now', term1
                              !WRITE(*,*) 'info for F', info
                              CALL ZGETRI(norb, term2, norb, 
     1                        ipiv, work, norb, info)
                              !WRITE(*,*) 'should be one', matmul(term1,part)
                              !WRITE(*,*) 'info for I', info
                              !WRITE(*,*) 'Haha'

                              DO l1 = 1, norb

                                  DO l2 = 1, norb

                                      DO l3 = 1, norb

                                          DO l4 = 1, norb

                                              term(l1,l2,l3,l4) = 
     1                                        term2(l3,l1)*term1(l2,l4) 

                                              chi_real(l1,l2,l3,l4,iq)=
     1                                        chi_real(l1,l2,l3,l4,iq)+
     1                                        REAL(term(l1,l2,l3,l4))

                                              chi_imag(l1,l2,l3,l4,iq)=
     1                                        chi_imag(l1,l2,l3,l4,iq)+
     1                                        AIMAG(term(l1,l2,l3,l4))

                                          ENDDO !DO l4

                                      ENDDO !DO l3
                                  
                                  ENDDO !DO l2

                              ENDDO !DO l1

                              !WRITE(*,*) 'Hahaha'                      

                          ENDDO !DO n
                      ENDIF !DO ikp
                  ENDDO !DO ikz
              ENDDO !DO iky
          ENDDO !DO ikx

C           ! the first block
C           DO d1=1,norb
C               l1 = d1
C               l2 = d1
C               DO d2=1,norb
C                   l3 = d2
C                   l4 = d2
C                   CHI_0(d1,d2,iq)=CMPLX(chi_real(l1,l2,l3,l4,iq),
C      1            chi_imag(l1,l2,l3,l4,iq))
C               ENDDO ! DO d2
C           ENDDO ! DO d1

C           DO c1=1,norb
C               index_r=1
C               DO c2=1,norb
C                   IF (c1 .NE. c2) THEN
C                       ! row index
C                       r_in = norb+(norb-1)*(c1-1)+index_r
C                       index_r = index_r + 1
C                       l1 = c1
C                       l2 = c2

C                       DO c3=1,norb
C                           index_c = 1
C                           DO c4=1,norb
C                               IF (c3 .NE. c4) THEN
C                                   !column index
C                                   c_in=norb+(norb-1)*(c3-1)+index_c
C                                   index_c = index_c+1
C                                   l3 = c3
C                                   l4 = c4
C                                   CHI_0(r_in,c_in,iq)=CMPLX(
C      1                            chi_real(l1,l2,l3,l4,iq),
C      1                            chi_imag(l1,l2,l3,l4,iq))
C                                ENDIF
C                           ENDDO ! DO C4
C                       ENDDO ! DO c3
C                   ENDIF
C               ENDDO ! DO c2
C           ENDDO ! DO c1

C           WRITE(*,*) 'Hello'

C           part = identity_S - matmul(U_S(:,:,iq),
C      1    (-T/numk)*CHI_0(:,:,iq))

C           CALL ZGETRF(norb*norb,norb*norb,part,norb*norb,ipiv,info)
C           CALL ZGETRI(norb*norb,part,norb*norb,ipiv,work,norb*norb,info)

C           CHI_S(:,:,iq)=matmul(part,(-T/numk)*CHI_0(:,:,iq))

          !WRITE(*,*) 'Flag up'
          !WRITE(*,*) chi_real(1,1,1,1,iq) 
      ENDDO !iq

      !DO i = 1, 60
          !GG_real(i) = GG_real(i)/numk
          !GG_complex(i) = GG_complex(i)/numk
      !ENDDO

      !CALL cpu_time(stop_time)

      !WRITE(*,*) stop_time - start_time, "seconds"

!####################################################################################################
      CALL MPI_REDUCE(chi_real,tot_chi_real,numq*norb*norb*norb*norb,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)
      CALL MPI_REDUCE(chi_imag,tot_chi_imag,numq*norb*norb*norb*norb,
     1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)
C       CALL MPI_REDUCE(CHI_S,tot_CHI_S,numq*norb*norb*norb*norb,
C      1MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_Comm_World,ierr)

      ! output
      IF (rank .EQ. 0) THEN
          DO i = 1, numq
              chi_real(:,:,:,:,i) = (-T/numk)*tot_chi_real(:,:,:,:,i)
              chi_imag(:,:,:,:,i) = (-T/numk)*tot_chi_imag(:,:,:,:,i)
          ENDDO

          DO iq = 1,numq
              ! build CHI_0
              ! the first block
              DO d1=1,norb
                  l1 = d1
                  l2 = d1
                  DO d2=1,norb
                      l3 = d2
                      l4 = d2
                      CHI_0(d1,d2,iq)=CMPLX(chi_real(l1,l2,l3,l4,iq),
     1                chi_imag(l1,l2,l3,l4,iq),8)
                  ENDDO ! DO d2
              ENDDO ! DO d1

              DO c1=1,norb
                  index_r=1
                  DO c2=1,norb
                      IF (c1 .NE. c2) THEN
                          ! row index
                          r_in = norb+(norb-1)*(c1-1)+index_r
                          index_r = index_r + 1
                          l1 = c1
                          l2 = c2

                          DO c3=1,norb
                              index_c = 1
                              DO c4=1,norb
                                  IF (c3 .NE. c4) THEN
                                      !column index
                                      c_in=norb+(norb-1)*(c3-1)+index_c
                                      index_c = index_c+1
                                      l3 = c3
                                      l4 = c4
                                      CHI_0(r_in,c_in,iq)=CMPLX(
     1                                chi_real(l1,l2,l3,l4,iq),
     1                                chi_imag(l1,l2,l3,l4,iq),8)
                                  ENDIF
                              ENDDO ! DO C4
                          ENDDO ! DO c3
                      ENDIF
                  ENDDO ! DO c2
              ENDDO ! DO c1

              part = identity_S - matmul(U_S,CHI_0(:,:,iq))

              CALL ZGETRF(norb*norb,norb*norb,part,norb*norb,ipivS,info)
              CALL ZGETRI(norb*norb,part,norb*norb,ipivS,workS,
     1        norb*norb,info)      

              CHI_S(:,:,iq)=matmul(part,CHI_0(:,:,iq))

              part = identity_S + matmul(U_S,CHI_0(:,:,iq))

              CALL ZGETRF(norb*norb,norb*norb,part,norb*norb,ipivS,info)
              CALL ZGETRI(norb*norb,part,norb*norb,ipivS,workS,
     1        norb*norb,info)      

              CHI_C(:,:,iq)=matmul(part,CHI_0(:,:,iq))

          ENDDO ! do iq

          DO i=1,numq
              WRITE(*,*) i, chi_real(1,1,1,1,i), chi_imag(1,1,1,1,i)
C               DO l1=1,norb*norb
C                   DO l2=1,norb*norb
C                       WRITE(*,*) l1,l2,CHI_0(l1,l2,i)
C                       WRITE(*,*) l1,l2,CHI_S(l1,l2,i)
C                   ENDDO
C               ENDDO
          ENDDO

      ENDIF !IF rank
      ! chi_real = (-T/numk)*tot_chi_real

      ! open(20, file='sigi.txt')

      ! close(20)
      DEALLOCATE(ommesh)

C       WRITE(*,*) 'Hehe'

      RETURN

      END SUBROUTINE     



