module const

  implicit none
  integer, parameter :: nx    = 2048       ! number of grid points
  integer, parameter :: nxgs  = 2         ! start point
  integer, parameter :: nxge  = nxgs+nx-1 ! end point
  integer, parameter :: np    = 1000     ! number of particles in each cell
  integer, parameter :: nsp   = 2         ! number of particle species
  integer, parameter :: nproc = 12         ! number of processors
  integer, parameter :: bc    = 0         ! boundary condition (periodic:0, reflective:-1)

end module
  
