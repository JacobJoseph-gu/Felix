
MODULE test_koch_mod

! equation 19
! A useful expansion of the exponential of the sum of two non-commuting matrices, one of which is diagonal
! By - Christoph T Koch and John C H Spence
! J. Phys. A: Math. Gen. 36 (2003) 803–816

! Dcoeff from equation 18

! bprime(0:u) are the unique B(l(0:q)), and there is a degeneracy d(k) for each bprime
! (bprime(k) unique -> d(k) = 0)
! bprime(k) used instead of bprime(l(k))

IMPLICIT NONE
CONTAINS

  ! S_n,m = e ** ( lambda ( A + B ) )
  ! where A, B matrices, B diagonal, lambda a scalar

  SUBROUTINE CalculateElementS ( lambda, A, Bmatrix, nnd, mnd, max_q, S )

    ! INPUTS
    INTEGER(4),INTENT(IN) :: nnd, mnd, max_q
    REAL(8),INTENT(IN) :: lambda, Bmatrix(:,:), A(:,:)
    ! OUTPUTS
    REAL(8),INTENT(OUT) :: S
    ! local variables
    REAL(8) :: Ccoeff, B(SIZE(A,1)),&
               sumproduct ! to store an iterative product
    INTEGER(4) :: ind,& ! a generic looping index
                  q, N, l(0:max_q)
    CHARACTER(50) :: formatting ! for writing terminal output

    IF(.NOT.(SIZE(A,1).EQ.SIZE(A,2).AND.SIZE(Bmatrix,1).EQ.SIZE(Bmatrix,2)&
          .AND.SIZE(A,1).EQ.SIZE(Bmatrix,1))) RETURN

    N = SIZE(A,1)

    ! PRINT A, B, lambda, 
    WRITE(*,'(a,i3,a,i3,a,F6.3)') " n = ", nnd, "     m = ", mnd, "    lamda = ", lambda 
    DO ind = 1,N
      IF(ind.EQ.1) THEN
        WRITE(formatting,'(a,i0,a,i0,a)') "(a,", N, "(F6.3,1x),a,", N, "(F6.3,1x),a)"
        WRITE(*,formatting) " A = [ ", A(:,ind), "]    B = [ ", Bmatrix(:,ind), "]"
      ELSE 
        WRITE(formatting,'(a,i0,a,i0,a)') "(a,", N, "(F6.3,1x),a,", N, "(F6.3,1x),a)"
        WRITE(*,formatting) "     [ ", A(:,ind), "]        [ ", Bmatrix(:,ind), "]"
      END IF 
    END DO
    WRITE(*,*) '-------------------------------------------------------------'

    ! use single dim N array for B instead of N x N diagonal matrix
    DO ind = 1,N
      B(ind) = Bmatrix(ind,ind)
    END DO
    S = 0

    ! + e ** ( lambda b_n ) * delta_n,m 
    IF(nnd.EQ.mnd) S = S + EXP(lambda*B(nnd))

    ! ------------------------------------
    ! summation over q
    ! ------------------------------------

    DO q = 1,max_q
      WRITE(*,'(a,i0)') 'q = ',q

      IF(q.EQ.1) THEN
        ! simply one term in summation a_n,m * Ccoeff
        l(0) = nnd
        l(1) = mnd
        CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff )
        S = S + A(nnd,mnd) * Ccoeff
      ELSEIF(q.EQ.2) THEN
        ! summation becomes a_n,l1 * a_l1,m * Ccoeff with l1 = 0 to l1 = N
        l(0) = nnd
        l(2) = mnd
        DO ind = 1,N
          l(1) = ind
          CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff )
          S = S + A(nnd,l(1)) * A(l(1),mnd) * Ccoeff
        END DO
      ELSE
        ! ------------------------------------
        ! q.GE.3
        ! ------------------------------------ 
        
        l = 1
        l(0) = nnd
        l(q) = mnd

        ! assume that can simply sum all of the final products together
        ! incrementing low indices l(1) first and l(q-1) last
  ! ---------------------------------------------------------------------------------
        DO WHILE (l(q-1).LT.N+1)

          !WRITE(*,'(a)') "--------------------------------------------------"
          !WRITE(formatting,'(a,i0,a)') "(a,", q-1, "(i4,1x))"
          !WRITE(*,formatting) 'l(1:q+1) = ', l(1:q+1)

          ! multiply a_n,l1 , a_l1,l2 , a_l2,l3, ... , a_lq-1,lq , a_lq,m
          ! A(l(0),l(1)) ... A(l(q-1),l(q))
          sumproduct = A(l(0),l(1))
          DO ind = 1,q-1
            sumproduct = sumproduct * A(l(ind),l(ind+1))
          END DO

          CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff ) 
          S = S + sumproduct * Ccoeff

          !WRITE(*,*) 'sumproduct'   
          !WRITE(*,'(F6.2)') Ccoeff
          !IF(ABS(Ccoeff).GT.1) WRITE(*,*) 'Ccoeff', Ccoeff, ' l(0:q) = ',l(0:q)
          !IF(q.EQ.7) WRITE(*,*) 'Ccoeff', Ccoeff,'||| q', q, ' l(0:q) = ',l(0:q)
          !WRITE(*,*) 'S',S

          ! iterate l_1 by 1 and check through each summation,
          ! if a summation has reached N, reset and increment summation above
          l(1) = l(1) + 1
          DO ind = 2,q-1
            IF(l(ind-1).EQ.N+1) THEN
              l(ind-1) = 1
              l(ind) = l(ind) + 1     
            ELSE
              EXIT
            END IF
          END DO

        END DO
  ! ---------------------------------------------------------------------------------
      END IF 
      WRITE(*,'(a,2x,F10.4,a)') 'S = ',S, ' -------------------------------------------------------------'  
    END DO

  END SUBROUTINE



  SUBROUTINE CalculateCcoeff ( B, lambda, l, q, Ccoeff )

    REAL(8),INTENT(OUT) :: Ccoeff
    REAL(8),INTENT(IN) :: B(:), lambda
    INTEGER(4),INTENT(IN) :: l(0:), q
    REAL(8) :: bprime(0:q), Dcoeff, iterationproduct
    INTEGER(4) :: d(0:q), u, jprime, k, r, ind
    CHARACTER(50) :: formatting ! for writing terminal output

    CALL GetUniqueSubset ( B, l, q, bprime, d, u )

    !WRITE(*,'(a)') "--------------------------------------------------"
    !WRITE(formatting,'(a,i0,a,i0,a)') "(a,", q+1, "(i4,1x)'|||'", q+1, "(F5.2,1x))"
    !WRITE(*,formatting) 'l(0:q), ( B(l(ind)), ind=1,q+1 ) ', l(0:q), ( B(l(ind)), ind=0,q )
    !WRITE(formatting,'(a,i0,a,i0,a)') '(a,', u+1, '(F5.2,1x)"|||",', u+1, '(i4,1x))'
    !WRITE(*,formatting) 'bprime(0:u), d(0:u)', bprime(0:u), d(0:u)
    !IF(q.EQ.3) THEN
    !  IF(ALL(l(0:q).EQ.[1,3,2,2])) THEN
    !    WRITE(formatting,'(a,i0,a,i0,a)') "(a,", q+1, "(i4,1x)'|||'", q+1, "(F5.2,1x))"
    !    WRITE(*,formatting) 'l(0:q), ( B(l(ind)), ind=1,q+1 ) ', l(0:q), ( B(l(ind)), ind=0,q )
    !    WRITE(formatting,'(a,i0,a,i0,a)') '(a,', u+1, '(F5.2,1x)"|||",', u+1, '(i4,1x))'
    !    WRITE(*,formatting) 'bprime(0:u), d(0:u)', bprime(0:u), d(0:u)
    !    WRITE(*,'(a)') '---------- SUM over k and jprime ----------------------------------- ' 
    !  END IF
    !END IF

    Ccoeff = 0
    DO k = 0,u
      DO jprime = 0,d(k)

        iterationproduct = 0
        ! summation from r = 0 to r = q - jprime - 1
        DO r = 0, q - jprime - 1

          !WRITE(*,*) lambda, bprime(k)
          !WRITE(*,'(E10.2)') (lambda * bprime(k))**REAL(r,kind(8)) / real(factorial(r),kind(8))
          
          iterationproduct = iterationproduct + (lambda * bprime(k))**REAL(r,kind(8)) / real(factorial(r),kind(8))
        END DO
        IF(q - jprime - 1.EQ.-1) Ccoeff = 0    
        iterationproduct = -iterationproduct

        !WRITE(*,*) 'iterationproduct only after sum over r', iterationproduct

        iterationproduct = iterationproduct + exp( lambda * bprime(k) )
        
        iterationproduct = iterationproduct * lambda**REAL(jprime,KIND(8)) / real(factorial(jprime),kind(8))

        CALL CalculateDcoeff ( B, l, q, bprime, k, d, u, jprime, Dcoeff )
        iterationproduct = iterationproduct * Dcoeff

        Ccoeff = Ccoeff + iterationproduct

        !IF(q.EQ.3) THEN
        !  IF(ALL(l(0:q).EQ.[1,3,2,2])) THEN
        !    WRITE(*,*) 'jprime = ',jprime, ' k = ',k, 'factorial & lambda', &
        !    lambda**REAL(jprime,KIND(8)) / real(factorial(jprime),kind(8))
        !    WRITE(*,*) 'l(0:q) ', l(0:q)        
        !    WRITE(*,*) 'Dcoeff', Dcoeff, 'k = ', k, ' jprime = ',jprime
        !    WRITE(*,'(a)') "--------------------------------------------------"
        !  END IF
        !END IF
        !WRITE(*,*) 'Dcoeff', Dcoeff, 'k = ', k, ' jprime = ',jprime
        !WRITE(*,*) ' exp( lambda * bprime(k) )', exp( lambda * bprime(k) )
        !WRITE(*,*) 'iterationproduct', iterationproduct
        !WRITE(*,*) 'k, jprime', k, jprime
        !WRITE(*,*) 'q - jprime - 1', q - jprime - 1
        !IF(ABS(Dcoeff).GT.1) WRITE(*,*) 'Dcoeff ', Dcoeff
        !IF(ABS(iterationproduct).GT.1) WRITE(*,*) 'iterationproduct ', iterationproduct

      END DO
    END DO

    !WRITE(*,*) 'Ccoeff', Ccoeff
    !IF(q.EQ.7) THEN
    !  IF(ALL(l(0:q).EQ.[1,2,2,2,2,1,1,1])) THEN
    !    WRITE(*,*) 'Ccoeff', Ccoeff
    !  END IF
    !END IF
    
  END SUBROUTINE CalculateCcoeff



  SUBROUTINE GetUniqueSubset ( B, l, q, bprime, d, u )

    REAL(8),INTENT(IN) :: B(:)
    REAL(8),INTENT(INOUT) :: bprime(0:)
    INTEGER(4),INTENT(INOUT) :: d(0:)
    INTEGER(4),INTENT(IN) :: l(0:), q
    INTEGER(4) :: ind, k,&
                  u ! number of unique elements minus 1
    LOGICAL :: IsUnique

    u = -1
    d = 0 ! degenearacy, b(l(k)) unique -> d(k) = 0
    DO ind = 0, q
      IsUnique = .TRUE.
      DO k = 0,u
        IF(B(l(ind)).EQ.bprime(k)) THEN
          IsUnique = .FALSE.
          d(k) = d(k) + 1
          EXIT
        END IF
      END DO
      IF(IsUnique) THEN
        bprime(u + 1) = B(l(ind))
        u = u + 1
      END IF         
    END DO

  END SUBROUTINE



  SUBROUTINE CalculateDcoeff ( B, l, q, bprime, k, d, u, jprime, Dcoeff )

    REAL(8),INTENT(IN) :: B(:), bprime(0:)
    INTEGER(4),INTENT(IN) :: l(0:), q, k, d(0:), u, jprime
    REAL(8),INTENT(OUT) :: Dcoeff
    INTEGER(4) :: ind, INoOfPermittedr, r_index
    INTEGER(4) :: r(d(k) - jprime), permitted_r_values(q+1), r_permitted_referance(d(k) - jprime)
    REAL(8) :: rsum(d(k) - jprime)
    LOGICAL :: productdefined, NotFinished
    CHARACTER(50) :: formatting ! for writing terminal output

    ! -----------------------------------------------------------
    ! sign (-1) ** and capitial pie product from r = 0 to r = u
    ! -----------------------------------------------------------
    
    Dcoeff = REAL( (-1) ** ( d(k) - jprime ) , KIND(8) )

    DO r_index = 0, u

      !WRITE(*,*) 'r_index, k', r_index, k

      IF(r_index.NE.k) THEN
        productdefined = .true.
        Dcoeff = Dcoeff * ( bprime(k) - bprime(r_index) ) ** (-REAL( d(r_index) + 1, KIND(8) ))

        !IF(q.EQ.7.AND.ALL(l(0:q).EQ.[1,2,2,2,2,1,1,1]).AND.k.EQ.0.AND.jprime.EQ.0) THEN
        !  WRITE(*,*) 'bprime(k) = ', bprime(k), ' bprime(r_index) = ',bprime(r_index)
        !  WRITE(*,*) 'fraction and d(r) power', ( bprime(k) - bprime(r_index) ) ** (-REAL( d(r_index) + 1, KIND(8) )) 
        !END iF

      END IF
    END DO

    IF(productdefined.EQV..FALSE.) THEN
      Dcoeff = 1
      RETURN
    END IF
    ! NB assumed that both this and embedded sum should have value = 1 for null case

    ! -----------------------------------------------------------
    ! large loop for r_1, ..., r_(d_k - j') summation & product
    ! -----------------------------------------------------------
  
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------        
    IF((d(k) - jprime).GT.0) THEN

      ! NB this always seems to lead to INoOfPermittedr.GT.0

      ! find permitted r values (this list will not contain duplicate r values)
      permitted_r_values = -1
      INoOfPermittedr = 0
      DO ind = 0,q

        !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
          !WRITE(*,*) 'B(l(ind))', B(l(ind)), '||| bprime(k)', bprime(k)
        !END iF
        !WRITE(*,*) 'B(l(ind)) , bprime(k)', B(l(ind)), '|||', bprime(k)

        IF(B(l(ind)).NE.bprime(k)) THEN
          INoOfPermittedr = INoOfPermittedr + 1
          permitted_r_values(INoOfPermittedr) = ind
        END IF
      END DO

      !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
        !WRITE(*,*) 'permitted_r_values', permitted_r_values
        !WRITE(*,*) 'd(k) - jprime', d(k) - jprime
      !END iF


      IF(INoOfPermittedr.GT.0) THEN

        r = permitted_r_values( 1 )
        r_permitted_referance =  1

        !WRITE(*,*) 'Iteration with INoOfPermittedr.GT.0'
        !WRITE(*,*) 'q, permitted_r_values', q, '|||', permitted_r_values
        !WRITE(*,*) 'r', r(:)
        !WRITE(formatting,'(a,i0,a,i0,a)') "(a,", q+1, "(i4,1x)','", q+1, "(F5.2,1x))"
        !WRITE(*,formatting) 'l(0:q), ( B(l(ind)), ind=0,q ) ', l(0:q), ( B(l(ind)), ind=0,q )
        !WRITE(formatting,'(a,i0,a,i0,a)') '(a,', u+1, '(F5.2,1x)",",', u+1, '(i4,1x))'
        !WRITE(*,formatting) 'bprime(0:u), d(0:u)', bprime(0:u), d(0:u)

        rsum = 0
        IF((d(k) - jprime).EQ.1) THEN
          DO ind = 1,INoOfPermittedr

            !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
              !WRITE(*,*) 'r', r, '||| all permitted ', permitted_r_values(1:INoOfPermittedr), '||| r ref ', r_permitted_referance
            !END iF
            !WRITE(*,*) 'r stuff', r, '|||', permitted_r_values(INoOfPermittedr), '|||', r_permitted_referance

            rsum(1) = rsum(1) + 1.0_8/(bprime(k) - B(l(r(1))))
            r(1) = permitted_r_values( ind+1 )
          END DO
        ELSE ! d(k) - jprime .GE. 2
  ! --------------------------------------------------------------------------------- 
          NotFinished = .TRUE.        
          DO WHILE (NotFinished)

            !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
            !  WRITE(*,*) 'r', r, '||| all permitted ', permitted_r_values(1:INoOfPermittedr), '||| r ref ', r_permitted_referance
            !END iF
            !WRITE(*,*) 'NotFinished',NotFinished
            !WRITE(*,*) 'r stuff', r, '|||', permitted_r_values(INoOfPermittedr), '|||', r_permitted_referance

            rsum(d(k) - jprime) = rsum(d(k) - jprime) + 1.0_8/(bprime(k) - B(l(r(d(k) - jprime))))

            ! moving outwards from innermost sum, check if max index of each sum has been reached
            DO ind = d(k) - jprime, 2, -1
              IF(r(ind).EQ.r(ind-1)) THEN
                ! contribute this summation to sum above
                rsum(ind-1) = rsum(ind-1) + 1.0_8/(bprime(k) - B(l(r(ind-1)))) * rsum(ind)          

                IF(ind.EQ.2) THEN
                  IF(r(1).LT.permitted_r_values(INoOfPermittedr)) THEN ! increment r(1)
                    r(1) = permitted_r_values( r_permitted_referance(1) + 1 )
                    r_permitted_referance(1) = r_permitted_referance(1) + 1
                  ELSE
                    NotFinished = .FALSE.
                  END IF
                END IF

                ! set this r to first permitted value
                rsum(ind) = 0
                r(ind) = permitted_r_values( 1 )
                r_permitted_referance(ind) = 1    
              ELSE ! r(ind).LT.r(ind-1)
                ! iterate sum to next permitted r value
                r(ind) = permitted_r_values( r_permitted_referance(ind) + 1 )
                r_permitted_referance(ind) = r_permitted_referance(ind) + 1      
                EXIT
              END IF
            END DO
          
          END DO
  ! ---------------------------------------------------------------------------------        
        END IF

        !IF(q.EQ.7.AND.ALL(l(0:q).EQ.[1,2,2,2,2,1,1,1]).AND.k.EQ.0.AND.jprime.EQ.0) THEN
        !  WRITE(*,*) 'rsum(1)', rsum(1)
        !END iF        
        !WRITE(*,*) rsum

        Dcoeff = Dcoeff * rsum(1)

      END IF

    END IF
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------         

  END SUBROUTINE



  RECURSIVE FUNCTION factorial ( n ) result ( f )
    INTEGER(4) :: f, n
    IF(n.EQ.1.OR.n.EQ.0) THEN
      f = 1
    ELSE
      f = n * factorial( n - 1_4 )
    END IF
  END FUNCTION

    

END MODULE
! to show in the meeting!
! test suite for range of values
! move inputs to function
! COMPLEX numbers or directly intensity
! impliment current debugging procedure
! impliment into felix
! HUGE PERFORMANCE IMPROVEMENTS
! max_q dynamical error handling
! impliment cleverer way to do permitted values using previous calculations, e.g. (q-1 - d(k)), u
! could impliment checking B diagonal, currently ignoring non-diagonal elements

! gfortran -o test_koch_mod test.f90 test_koch_mod.f90 -fbounds-check
! doing clever permutation stuff for non-unique multiplication
! google whether this has been programmed before
! affect of real byte size on performance
! internal functions have access to external scope but can overwrite
! variable scope for clarity?
! pure/elemental procedures
! use pointers to referance B elements?
! B diagonal so store as single dim array
