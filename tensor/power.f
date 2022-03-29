      PROGRAM POWER
      ! IMPLICIT NONE
      REAL*8, parameter :: pi = 3.14159265358979

      REAL*8, ALLOCATABLE :: xr(:, :), xi(:, :)      
      COMPLEX*16, ALLOCATABLE :: A(:, :), x_k(:, :), x_kp1(:, :)
      COMPLEX*16 :: leigenvalue
      INTEGER :: dim, i, num_iter

      REAL*8 :: start_time, stop_time

      ! for lapack
      COMPLEX*16, ALLOCATABLE :: eigenvalues(:)
      REAL*8, ALLOCATABLE :: eigenvalues_real(:)
      INTEGER :: LDVL, LDVR, INFO, LWORK
      INTEGER :: location
      COMPLEX*16, ALLOCATABLE :: VL(:, :), VR(:, :), WORK(:)
      REAL*8, ALLOCATABLE :: RWORK(:) 

      dim = 1000000
      ALLOCATE(A(dim, dim))
C       ALLOCATE(xr(dim, 1))
C       ALLOCATE(xi(dim, 1))
C       ALLOCATE(x_k(dim, 1))
C       ALLOCATE(x_kp1(dim, 1))
      
      WRITE(*,*) 'Hello'

      LDVL = dim
      LDVR = dim
      ALLOCATE(eigenvalues(dim))
      ALLOCATE(eigenvalues_real(dim))
      ALLOCATE(VL(dim, dim))
      ALLOCATE(VR(dim, dim))

      LWORK = 2*dim
      ALLOCATE(WORK(LWORK))
      ALLOCATE(RWORK(LWORK))

      WRITE(*,*) 'Hey'

C       num_iter = 500

      A = 0.0
C       CALL RANDOM_NUMBER(xr)
C       CALL RANDOM_NUMBER(xi)
C       x_k = CMPLX(xr, xi)
C       x_kp1 = 0.0

      WRITE(*,*) 'Ha'

      ! construct A
      DO i = 1, dim-1
         A(i, i) = 1
      ENDDO
      A(dim, dim) = 100000

      !CALL cpu_time(start_time)

C       ! power method
C       DO i=1, num_iter

C          x_kp1 = matmul(A, x_k)
C          x_k = x_kp1 / (NORM2(REAL(x_kp1))+NORM2(AIMAG(x_kp1)))

C       ENDDO

C       leigenvalue = dot_product(RESHAPE(matmul(transpose(x_k), A), 
C      1(/dim/)), RESHAPE(x_k, (/dim/)))

C       WRITE(*, *) 'The largest eigenvalue:', leigenvalue
C       CALL cpu_time(stop_time)
C       WRITE(*,*) 'Total time is ', stop_time - start_time, "seconds"

C       WRITE(*,*) 'Haha'

      CALL cpu_time(start_time)

      ! regular lapack subroutine
      CALL ZGEEV('N', 'N', dim, A, dim, eigenvalues, VL, LDVL, VR, 
     1LDVR, WORK, LWORK, RWORK, INFO)

      WRITE(*,*) 'Hahaha'

      eigenvalues_real = ABS(eigenvalues)
      location = MAXLOC(eigenvalues_real, 1)
      WRITE(*,*) 'location:', location

      leigenvalue = eigenvalues(location)

      WRITE(*, *) 'The largest eigenvalue:', leigenvalue

      CALL cpu_time(stop_time)
      WRITE(*,*) 'Total time is ', stop_time - start_time, "seconds"

      DEALLOCATE(A)
C       DEALLOCATE(xr)
C       DEALLOCATE(xi)
C       DEALLOCATE(x_k)
C       DEALLOCATE(x_kp1)

      DEALLOCATE(eigenvalues)
      DEALLOCATE(eigenvalues_real)
      DEALLOCATE(VL)
      DEALLOCATE(VR)
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)

      END PROGRAM POWER
