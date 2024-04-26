      module PCG

      integer :: N2
      integer :: NLUmax, NLU, METHOD
      integer :: NPLU

      real(kind=8) :: EPSICCG

      real(kind=8), dimension(:), allocatable :: D, PHI, BFORCE
      real(kind=8), dimension(:), allocatable :: AMAT

      integer, dimension(:), allocatable :: INLU
      integer, dimension(:), allocatable :: indexLU, itemLU

      end module PCG
