      PROGRAM DMFT_chi_static_sp
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
      COMPLEX*16,ALLOCATABLE :: sigi_up(:,:)! self-energy, spin up
      COMPLEX*16,ALLOCATABLE :: sigi_down(:,:)! self-energy, spin down
      REAL*8,ALLOCATABLE :: chi_real(:) ! static irreducible susceptibility
      REAL*8 :: distance  ! distance of q-points

      REAL*8  :: numesh ! bosonic Matsubara frequency
      INTEGER :: i,j,k    ! dummy index
      LOGICAL :: iffile

      ! for MPI
      INTEGER rank,total,ierr,errorcode

      WRITE(*,*) 'DMFT_chi_static_sp starts:'

      INQUIRE(FILE = 'DMFT_chi_static_sp.in', EXIST = iffile)
      IF (iffile .EQV. .false.) THEN
          WRITE(*,*) 'input is needed for DMFT_chi_static_sp!'
          STOP
      ELSE
          OPEN(65, FILE = "DMFT_chi_static_sp.in")

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
          ALLOCATE(sigi_up(ncor_orb,2*Nc))
          READ(65,*) ((sigi_up(i,j),j=1,2*Nc),i=1,ncor_orb)
          ALLOCATE(sigi_down(ncor_orb,2*Nc))
          READ(65,*) ((sigi_down(i,j),j=1,2*Nc),i=1,ncor_orb)
          CLOSE(65)

      ENDIF !IF iffile

      ALLOCATE(chi_real(numq))

      ! launch MPI
      CALL MPI_INIT(ierr)

      IF (ierr .NE. MPI_SUCCESS) THEN
          WRITE(*,*), 'Error starting MPI program. Terminating.'
          CALL MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      ENDIF !IF ierr

      CALL MPI_Comm_size(MPI_Comm_World, total, ierr)
      CALL MPI_Comm_rank(MPI_Comm_World, rank, ierr)

      CALL do_chi_static_sp(numq,q,T,Nc,mu,nk,nr,tran,ham,
     1sigi_up,sigi_down,chi_real,rank,total,norb,ncor_orb)

      IF (rank .eq. 0) THEN

          OPEN(65, FILE = "DMFT_chi_static_sp.out")

          distance = 0.0

          DO iq = 1, numq

              WRITE(65,'(5F24.16)') q(iq,1),q(iq,2),q(iq,3),
     1        distance,chi_real(iq)
              !WRITE(*,*) q(iq,1),q(iq,2),q(iq,3),distance,chi_real(iq)

              IF (iq .lt. numq) THEN
                  distance = distance + sqrt((q(iq+1,1) - q(iq,1))**2 
     1            + (q(iq+1,2) - q(iq,2))**2 + (q(iq+1,3) - q(iq,3))**2) 
              ENDIF

          ENDDO !iq

          CLOSE(65)

      ENDIF !IF rank
  
      DEALLOCATE(chi_real)   
  
      CALL MPI_FINALIZE(ierr)

      END PROGRAM DMFT_chi_static_sp
      