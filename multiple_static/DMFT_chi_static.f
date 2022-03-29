      PROGRAM DMFT_chi_static
      !IMPLICIT NONE
      INCLUDE 'mpif.h'
      REAL*8, PARAMETER :: pi = 3.14159265358979

      INTEGER :: numq   ! number of q points
      INTEGER :: iq     ! index of q points
      REAL*8,ALLOCATABLE  :: q(:,:)      ! q point
      REAL*8  :: T      ! temperature
      INTEGER :: Nc     ! the cutoff for fermionic Matsubara frequency
      REAL*8  :: mu     ! chemical potential
      INTEGER :: nk(3) ! number of k-point
      INTEGER :: nr                    ! number of "r" points

      INTEGER :: norb,ncor_orb

      INTEGER,ALLOCATABLE :: tran(:,:) ! R
      COMPLEX*16,ALLOCATABLE :: ham(:,:,:) ! H(R)
      COMPLEX*16,ALLOCATABLE :: sigi(:,:)! self-energy
      COMPLEX*16,ALLOCATABLE :: U_S(:,:) ! to calculate chi_S
      COMPLEX*16,ALLOCATABLE :: CHI_S(:,:,:)
      COMPLEX*16,ALLOCATABLE :: CHI_C(:,:,:)
      REAL*8,ALLOCATABLE :: chi_real(:,:,:,:,:) ! static irreducible susceptibility
      REAL*8 :: distance  ! distance of q-points
      REAL*8 :: trace,trace_S,trace_C

      !REAL*8  :: numesh ! bosonic Matsubara frequency
      INTEGER :: i,j,k    ! dummy index
      LOGICAL :: iffile

      ! for MPI
      INTEGER rank,total,ierr,errorcode

      WRITE(*,*) 'DMFT_chi_static starts:'

      INQUIRE(FILE = 'DMFT_chi_static.in', EXIST = iffile)
      IF (iffile .EQV. .false.) THEN
          WRITE(*,*) 'input is needed for DMFT_chi_static!'
          STOP
      ELSE
          OPEN(65, FILE = "DMFT_chi_static.in")

          READ(65,*) numq
          !WRITE(*,*) 'numq', numq

          ALLOCATE(q(numq,3))
          READ(65,*) ((q(i,j),j=1,3),i=1,numq)
          !WRITE(*,*) 'q', ((q(i,j),j=1,3),i=1,numq)

          READ(65,*) T
          !WRITE(*,*) 'T', T
          READ(65,*) Nc
          !WRITE(*,*) 'Nc', Nc
          READ(65,*) mu
          !WRITE(*,*) 'mu', mu
          READ(65,*) (nk(i),i=1,3)
          !WRITE(*,*) 'nk', (nk(i),i=1,3)

          READ(65,*) nr
          !WRITE(*,*) 'nr', nr
          READ(65,*) norb
          !WRITE(*,*) 'norb', norb
          READ(65,*) ncor_orb !sigma's dimension

          ALLOCATE(tran(nr,3))
          READ(65,*) ((tran(i,j),j=1,3),i=1,nr)
          !WRITE(*,*) 'tran', ((tran(i,j),j=1,3),i=1,nr)

          ALLOCATE(ham(nr,norb,norb))
          READ(65,*) (((ham(i,j,k),k=1,norb),j=1,norb),i=1,nr)
          !WRITE(*,*) 'ham', (((ham(i,j,k),k=1,norb),j=1,norb),i=1,nr)
          ALLOCATE(sigi(ncor_orb,2*Nc))
          READ(65,*) ((sigi(i,j),j=1,2*Nc),i=1,ncor_orb)

          ALLOCATE(U_S(norb*norb,norb*norb))
          READ(65,*) ((U_S(i,j),j=1,norb*norb),i=1,norb*norb)
C           U_S = 0.0
C           DO i=1,norb*norb
C               U_S(i,i) = 1.0
C           ENDDO

          CLOSE(65)

      ENDIF !IF iffile

      ALLOCATE(chi_real(norb,norb,norb,norb,numq))
      ALLOCATE(CHI_S(norb*norb,norb*norb,numq))
      ALLOCATE(CHI_C(norb*norb,norb*norb,numq))

      ! launch MPI
      CALL MPI_INIT(ierr)

      IF (ierr .NE. MPI_SUCCESS) THEN
          WRITE(*,*), 'Error starting MPI program. Terminating.'
          CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      ENDIF !IF ierr

      CALL MPI_Comm_size(MPI_Comm_World, total, ierr)
      CALL MPI_Comm_rank(MPI_Comm_World, rank, ierr)

      CALL do_chi_static(numq,q,T,Nc,mu,nk,nr,tran,ham,
     1sigi,U_S,CHI_S,CHI_C,chi_real,rank,total,norb,ncor_orb)

C       WRITE(*,*) 'Hi'      

      IF (rank .eq. 0) THEN

C           WRITE(*,*) 'Hello'

          OPEN(65, FILE = "DMFT_chi_static.out")

          distance = 0.0

          DO iq = 1, numq

C               WRITE(*,*) 'Hey'
              trace = 0.0
              trace_S = 0.0
              trace_C = 0.0

              DO i=1,norb
                  trace = trace + chi_real(i,i,i,i,iq)
                  trace_S = trace_S + REAL(CHI_S(i,i,iq))
                  trace_C = trace_C + REAL(CHI_C(i,i,iq))
              ENDDO

              WRITE(65,'(5F24.16)') q(iq,1),q(iq,2),q(iq,3),
     1        distance,trace
              WRITE(65,'(4F24.16)') chi_real(1,1,1,1,iq),
     1        chi_real(2,2,2,2,iq),chi_real(3,3,3,3,iq),
     1        chi_real(4,4,4,4,iq)
              WRITE(65,'(5F24.16)') REAL(CHI_S(1,1,iq)),
     1        REAL(CHI_S(2,2,iq)),REAL(CHI_S(3,3,iq)),
     1        REAL(CHI_S(4,4,iq)),trace_S 
              WRITE(65,'(5F24.16)') REAL(CHI_C(1,1,iq)),
     1        REAL(CHI_C(2,2,iq)),REAL(CHI_C(3,3,iq)),
     1        REAL(CHI_C(4,4,iq)),trace_C 
              !WRITE(*,*) q(iq,1),q(iq,2),q(iq,3),distance,chi_real(iq)

C               WRITE(*,*) 'Ha'

              IF (iq .lt. numq) THEN
                  distance = distance + sqrt((q(iq+1,1) - q(iq,1))**2 
     1            + (q(iq+1,2)-q(iq,2))**2 + (q(iq+1,3) - q(iq,3))**2) 
              ENDIF

          ENDDO !iq

          CLOSE(65)

      ENDIF !IF rank
      
      DEALLOCATE(q)
      DEALLOCATE(tran)
      DEALLOCATE(ham)
      DEALLOCATE(sigi)
      DEALLOCATE(chi_real)
      DEALLOCATE(CHI_S)
      DEALLOCATE(CHI_C)
      DEALLOCATE(U_S)   
  
      CALL MPI_FINALIZE(ierr)

      END PROGRAM DMFT_chi_static


