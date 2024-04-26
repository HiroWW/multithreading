!C
!C***
!C*** POI_GEN
!C***
!C
!C    generate COEF. MATRIX for POISSON equations
!C    
      subroutine POI_GEN

      use STRUCT
      use PCG

      implicit REAL*8 (A-H,O-Z)

!C
!C-- INIT.
      nn = ICELTOT

      NLU= 6

      allocate (BFORCE(nn), D(nn), PHI(nn))
      allocate (INLU(nn))

      PHI   = 0.d0
      BFORCE= 0.d0
        D   = 0.d0

      INLU= 0

!C
!C-- INTERIOR & NEUMANN boundary's
      do icel= 1, ICELTOT
        icN1= NEIBcell(icel,1)
        icN2= NEIBcell(icel,2)
        icN3= NEIBcell(icel,3)
        icN4= NEIBcell(icel,4)
        icN5= NEIBcell(icel,5)
        icN6= NEIBcell(icel,6)

        if (icN5.ne.0) then
          INLU(icel)= INLU(icel) + 1
        endif

        if (icN3.ne.0) then
          INLU(icel)= INLU(icel) + 1           
        endif

        if (icN1.ne.0) then
          INLU(icel)= INLU(icel) + 1           
        endif

        if (icN2.ne.0) then
          INLU(icel)= INLU(icel) + 1           
        endif

        if (icN4.ne.0) then
          INLU(icel)= INLU(icel) + 1           
        endif

        if (icN6.ne.0) then
          INLU(icel)= INLU(icel) + 1           
        endif

      enddo

!C
!C-- 1D array
      allocate (indexLU(0:nn))
      indexLU= 0
     
      do icel= 1, ICELTOT
        indexLU(icel)= INLU(icel)
      enddo

      do icel= 1, ICELTOT
        indexLU(icel)= indexLU(icel) + indexLU(icel-1)
      enddo

      NPLU= indexLU(ICELTOT)

      allocate (itemLU(NPLU), AMAT(NPLU))

      itemLU= 0
      AMAT= 0.d0
!C===

!C
!C +-----------------------------------+
!C | INTERIOR & NEUMANN BOUNDARY CELLs |
!C +-----------------------------------+
!C===
      icouG= 0
      do icel= 1, ICELTOT
        icN1= NEIBcell(icel,1)
        icN2= NEIBcell(icel,2)
        icN3= NEIBcell(icel,3)
        icN4= NEIBcell(icel,4)
        icN5= NEIBcell(icel,5)
        icN6= NEIBcell(icel,6)


        VOL0= VOLCEL(icel)

        icou= 0
        if (icN5.ne.0) then
          coef   =RDZ * ZAREA
          D(icel)= D(icel) - coef
                   icou= icou + 1
                   k   = icou + indexLU(icel-1)
          itemLU(k)= icN5
            AMAT(k)= coef
        endif

        if (icN3.ne.0) then
          coef   = RDY * YAREA
          D(icel)= D(icel) - coef
                   icou= icou + 1
                   k   = icou + indexLU(icel-1)
          itemLU(k)= icN3
            AMAT(k)= coef
        endif

        if (icN1.ne.0) then
          coef   = RDX * XAREA
          D(icel)= D(icel) - coef
                   icou= icou + 1        
                   k   = icou + indexLU(icel-1)
          itemLU(k)= icN1
            AMAT(k)= coef
        endif

        if (icN2.ne.0) then
          coef   = RDX * XAREA
          D(icel)= D(icel) - coef
                   icou= icou + 1
                   k   = icou + indexLU(icel-1)
          itemLU(k)= icN2
            AMAT(k)= coef
        endif

        if (icN4.ne.0) then
          coef   = RDY * YAREA
          D(icel)= D(icel) - coef
                   icou= icou + 1
                   k   = icou + indexLU(icel-1)
          itemLU(k)= icN4
            AMAT(k)= coef
        endif

        if (icN6.ne.0) then
          coef   = RDZ * ZAREA
          D(icel)= D(icel) - coef
                   icou= icou + 1
                   k   = icou + indexLU(icel-1)
          itemLU(k)= icN6
            AMAT(k)= coef
        endif

        ii= XYZ(icel,1)
        jj= XYZ(icel,2)
        kk= XYZ(icel,3)

        BFORCE(icel)= -dfloat(ii+jj+kk) * VOL0

      enddo
!C===

!C
!C +--------------------------+
!C | DIRICHLET BOUNDARY CELLs |
!C +--------------------------+
!C   TOP SURFACE
!C===
      do ib= 1, ZmaxCELtot
        icel= ZmaxCEL(ib)
        coef= 2.d0 * RDZ * ZAREA
        D(icel)= D(icel) - coef
      enddo
!C===

      return
      end
