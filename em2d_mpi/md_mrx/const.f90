module const

  implicit none

!!************************ NUMERICAL CONSTANTS ***********************************!!
  integer, parameter :: nx   = 320+1     ! NUMBER OF GRID POINTS IN X
  integer, parameter :: ny   = 640       ! NUMBER OF GRID POINTS IN Y
  integer, parameter :: nxgs = 2         ! START POINT IN X
  integer, parameter :: nxge = nxgs+nx-1 ! END POINT
  integer, parameter :: nygs = 2         ! START POINT IN Y
  integer, parameter :: nyge = nygs+ny-1 ! END POINT
  integer, parameter :: np   = 200*nx    ! CUMULATIVE NUMBER OF PARTICLES IN EACH Y POSITION
  integer, parameter :: nsfo = 1         ! SHAPE FUNCTION ORDER (0:NGP, 1:CIC, 2:SPLINE)
  integer, parameter :: nsp  = 2         ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: bc   = -1        ! BOUNDARY CONDITION (PERIODIC:0, REFLECTIVE:-1)
  integer, parameter :: nproc = 12       ! NUMBER OF PROCESSORS

!! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  integer, parameter :: itmax  = 20000	!NUMBER OF ITERATION
  integer            :: it0    = 0	!0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 500	!INTERVAL FOR PARTICLES & FIELDS STORAGE
  integer, parameter :: intvl2 = 50	!INTERVAL FOR ENERGY CALC.
  character(len=128) :: dir    = './dat/'	!DIRECTORY FOR OUTPUT
  character(len=128) :: fname_param  = 'init_param.dat'	!FILE NAME OF INIT CONDITIONS
  character(len=128) :: fname_energy = 'energy.dat'	!FILE NAME OF PARTICLE DATA

!! OTHER CONSTANTS
  real(8), parameter :: gfac   = 0.501D0 !IMPLICITNESS FACTOR > 0.5
  real(8), parameter :: cfl    = 0.5D0   !CFL CONDITION FOR LIGHT WAVE
  real(8), parameter :: delx   = 1.0D0   !CELL WIDTH
  real(8), parameter :: rdbl   = 1.0D0   !DEBYE LENGTH / CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*atan(1.0D0)

!!************************ PHYSICAL CONSTANTS ***********************************!!
!!       c : SPEED OF LIGHT
!!      mr : ION-TO-ELECTRON MASS RATIO
!!   alpha : wpe/wge = c/vth_e * sqrt(beta_e)
!!    beta : ION PLASMA BETA
!!   rtemp : Te/Ti
!! ncs,ngb : HARRIS SHEET DENSITY AND THE BACKGROUND DENSITY
!!     lcs : CURRENT SHEET THICKNESS IN UNIT OF ION SKIN DEPTH
  real(8), parameter :: c      = 1.0D0
  real(8), parameter :: mr     = 16.0D0
  real(8), parameter :: alpha  = 2.0D0, beta = 0.5D0, rtemp=0.25D0
  integer, parameter :: ncs    = 250, nbg  =  50
  real(8)            :: lcs    = 0.5D0

end module
