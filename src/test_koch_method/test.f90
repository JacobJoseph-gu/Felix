program main

  use test_koch_mod
  implicit none

  ! INPUTS
  INTEGER(4) :: nnd, mnd, max_q
  REAL(8),ALLOCATABLE :: lambda, Bmatrix(:,:), A(:,:)
  ! OUTPUTS
  REAL(8) :: S

  ! local
  INTEGER(4) :: N

  N = 3
  ALLOCATE( A(N,N) )
  ALLOCATE( Bmatrix(N,N) )

  lambda = 1
  A(:,1) =        [1,1,1]
  A(:,2) =        [1,1,1]
  A(:,3) =        [1,1,1]
  Bmatrix(:,1) =  [1,0,0]
  Bmatrix(:,2) =  [0,2,0]
  Bmatrix(:,3) =  [0,0,3]
  nnd = 1
  mnd = 2
  max_q = 10

  call CalculateElementS( lambda, A, Bmatrix, nnd, mnd, max_q, S )

  DEALLOCATE ( A, Bmatrix )

  N = 4
  ALLOCATE( A(N,N) )
  ALLOCATE( Bmatrix(N,N) )

  lambda = 1.1
  A(:,1) =        [0.01,0.01,0.01,0.02]
  A(:,2) =        [0.01,0.03,0.01,0.01]
  A(:,3) =        [0.01,0.01,0.02,0.01]
  A(:,4) =        [0.01,0.01,0.02,0.01]
  Bmatrix(:,1) =  [1.0,0.0,0.0,0.0]
  Bmatrix(:,2) =  [0.0,2.0,0.0,0.0]
  Bmatrix(:,3) =  [0.0,0.0,1.0,0.0]
  Bmatrix(:,4) =  [0.0,0.0,0.0,1.3]
  nnd = 2
  mnd = 2
  max_q = 7

  call CalculateElementS( lambda, A, Bmatrix, nnd, mnd, max_q, S )

  DEALLOCATE ( A, Bmatrix )

end program
