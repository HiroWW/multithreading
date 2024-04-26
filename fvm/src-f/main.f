      program MAIN

      use STRUCT
      use PCG
      use solver_PCG

      implicit REAL*8 (A-H,O-Z)

!C     
!C-- INIT.
      call INPUT
      call POINTER_INIT
      call BOUNDARY_CELL
      call CELL_METRICS
      call POI_GEN

!C 
!C-- MAIN SOLVER
      PHI=  0.d0

      ISET= 0

        call solve_PCG                                                  &
     &      ( ICELTOT, NPLU, indexLU, itemLU,  D, BFORCE,  PHI, AMAT,   &
     &     EPSICCG, ITR, IER)

      write (*,'(//, a,i10,1pe16.6)') '##ANSWER', ICELTOT, PHI(ICELTOT)        

      call OUTUCD

      stop
      end
