!C*** 
!C*** module solver_PCG
!C***
!
      module solver_PCG
      contains
!C
!C*** solve_PCG
!C
      subroutine solve_PCG                                              &
     &         ( N, NPLU, indexLU, itemLU, D, B, X, AMAT, EPS, ITR, IER)

      implicit REAL*8 (A-H,O-Z)

      real(kind=8), dimension(N)   :: D
      real(kind=8), dimension(N)   :: B
      real(kind=8), dimension(N)   :: X
      real(kind=8), dimension(NPLU):: AMAT

      integer, dimension(0:N) :: indexLU
      integer, dimension(NPLU)::  itemLU

      real(kind=8), dimension(:,:), allocatable :: W

      integer, parameter ::  R= 1
      integer, parameter ::  Z= 2
      integer, parameter ::  Q= 2
      integer, parameter ::  P= 3
      integer, parameter :: DD= 4

!C
!C +------+
!C | INIT |
!C +------+
!C===
      allocate (W(N,4))


      do i= 1, N
        X(i)  = 0.d0
        W(i,1)= 0.0D0
        W(i,2)= 0.0D0
        W(i,3)= 0.0D0
        W(i,4)= 0.0D0
      enddo

      do i= 1, N
        W(i,DD)= 1.d0/D(i)
      enddo
!C===

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===
      do i= 1, N
        VAL= D(i)*X(i)
        do k= indexLU(i-1)+1, indexLU(i)
          VAL= VAL + AMAT(k)*X(itemLU(k))
        enddo 
        W(i,R)= B(i) - VAL
      enddo

      BNRM2= 0.0D0
      do i= 1, N
        BNRM2 = BNRM2 + B(i)  **2
      enddo
!C===

!C
!C***************************************************************  ITERATION
      ITR= N

      do L= 1, ITR
!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===      
      do i= 1, N
        W(i,Z)= W(i,R)*W(i,DD)
      enddo
!C===

!C
!C +-------------+
!C | RHO= {r}{z} |
!C +-------------+
!C===
      RHO= 0.d0
      do i= 1, N
        RHO= RHO + W(i,R)*W(i,Z)   
      enddo          
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( L.eq.1 ) then

        do i= 1, N
          W(i,P)= W(i,Z)
        enddo
       else
         BETA= RHO / RHO1
         do i= 1, N
           W(i,P)= W(i,Z) + BETA*W(i,P)
         enddo
      endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===        
      do i= 1, N
        VAL= D(i)*W(i,P)
        do k= indexLU(i-1)+1, indexLU(i)
          VAL= VAL + AMAT(k)*W(itemLU(k),P)
        enddo 
        W(i,Q)= VAL
      enddo
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C1= 0.d0
      do i= 1, N
        C1= C1 + W(i,P)*W(i,Q)
      enddo

      ALPHA= RHO / C1
!C===


!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===
      do i= 1, N
        X(i)  = X(i)   + ALPHA * W(i,P)
        W(i,R)= W(i,R) - ALPHA * W(i,Q)
      enddo

      DNRM2= 0.d0
      do i= 1, N
        DNRM2= DNRM2 + W(i,R)**2
      enddo
!C===

      ERR = dsqrt(DNRM2/BNRM2)
      if (mod(L,100).eq.1) then
        write (*,'(i5,2(1pe16.6))') L, ERR
      endif

        if (ERR .lt. EPS) then
          IER = 0
          goto 900
         else
          RHO1 = RHO
        endif

      enddo
      IER = 1

  900 continue

      write (*,'(i5,2(1pe16.6))') L, ERR
      ITR= L
!      EPS= ERR

      deallocate (W)

      return

      end subroutine  solve_PCG
      end module     solver_PCG



