module const

  implicit none
  integer, parameter :: nx    = 128       ! number of grid points in x
  integer, parameter :: ny    = 128       ! number of grid points in y
  integer, parameter :: nxgs  = 2         ! start point in x
  integer, parameter :: nxge  = nxgs+nx-1 ! end point
  integer, parameter :: nygs  = 2         ! start point in y
  integer, parameter :: nyge  = nygs+ny-1 ! end point
  integer, parameter :: np    = 150*nx    ! number of particles in each cell
  integer, parameter :: nsp   = 2         ! number of particle species
  integer, parameter :: nproc = 8         ! number of processors
  integer, parameter :: bc    = 0         ! boundary condition in x (0:periodic, -1:reflective)

end module
  
