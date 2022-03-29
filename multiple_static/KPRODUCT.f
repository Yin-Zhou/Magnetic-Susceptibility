      SUBROUTINE KPRODUCT(A,B,AB)
      ! this subroutine calculates the direct product of A and B
      ! A and B are both two-dimensional arrays

      !IMPLICIT NONE

      COMPLEX*16 :: A(:, :), B(:, :)
      COMPLEX*16, ALLOCATABLE :: AB(:, :)

      INTEGER :: RA, RB, CA, CB
      INTEGER :: i, j, R

      RA = UBOUND(A, DIM=1)
      CA = UBOUND(A, DIM=2)
      RB = UBOUND(B, DIM=1)
      CB = UBOUND(B, DIM=2)

      IF (ALLOCATED(AB)) THEN
          DEALLOCATE(AB)
      ENDIF

      ALLOCATE(AB(RA*RB, CA*CB))

      R = 0

      DO i = 1, RA

          C = 0

          DO j = 1, CA

              AB((R+1):(R+RB), (C+1):(C+CB)) = A(i,j)*B
              C = C + CB

          ENDDO

          R = R + RB

      ENDDO

      RETURN

      END

