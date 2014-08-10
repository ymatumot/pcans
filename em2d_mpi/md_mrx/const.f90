module const

  implicit none
!  integer, parameter :: nx    = 406+1     ! number of grid points in x
!  integer, parameter :: ny    = 812       ! number of grid points in y
  integer, parameter :: nx    = 320+1     ! number of grid points in x
  integer, parameter :: ny    = 640       ! number of grid points in y
  integer, parameter :: nxgs  = 2         ! start point in x
  integer, parameter :: nxge  = nxgs+nx-1 ! end point
  integer, parameter :: nygs  = 2         ! start point in y
  integer, parameter :: nyge  = nygs+ny-1 ! end point
  integer, parameter :: np    = 100*nx    ! number of particles in each cell
  integer, parameter :: nsp   = 2         ! number of particle species
  integer, parameter :: nsfo  = 2         ! shape function order (0:NGP, 1:CIC, 2:Spline)
  integer, parameter :: nproc = 64        ! number of processors
  integer, parameter :: bc    = -1        ! boundary condition in x (0:periodic, -1:reflective)

end module
