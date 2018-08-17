module const

  implicit none

!!************************ NUMERICAL CONSTANTS ***********************************!!
  integer, parameter :: nx  = 2048  ! NUMBER OF GRID POINTS
  integer, parameter :: np  = 1024  ! NUMBER OF PARTICLES IN EACH CELL
  integer, parameter :: nsp = 2     ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: bc  = -1    ! BOUNDARY CONDITION (PERIODIC:0, REFLECTIVE:-1)

!! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  integer, parameter :: itmax  = 40000 !NUMBER OF ITERATION
  integer            :: it0    = 0     !0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 5000  !INTERVAL FOR PARTICLES & FIELDS STORAGE
  integer, parameter :: intvl3 = 1000  !INTERVAL FOR RECORDING MOMENT & FIELDS DATA	
  integer, parameter :: intvl4 = 1000  !INTERVAL FOR RECORDING PSD DATA
  character(len=128) :: dir    = './dat/' !DIRECTORY FOR OUTPUT
  character(len=128) :: dir_mom= './mom/' !DIRECTORY FOR MOMENT DATA
  character(len=128) :: dir_psd= './psd/' !DIRECTORY FOR MOMENT DATA
  character(len=128) :: file9  = 'init_param.dat' !FILE NAME OF INIT CONDITIONS
  character(len=128) :: file10 = 'file10.dat' !FILE NAME OF PARTICLE DATA

!! OTHER CONSTANTS
  real(8), parameter :: gfac   = 0.501D0 !IMPLICITNESS FACTOR > 0.5
  real(8), parameter :: cfl    = 0.5D0   !CFL CONDITION FOR LIGHT WAVE
  real(8), parameter :: rdbl   = 1.0D0   !DEBYE LENGTH / CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*atan(1.0D0)

!!************************ PHYSICAL CONSTANTS ***********************************!!
!!      n0 : NUMBER OF PARTICLES/CELL IN THE UPSTREAM REGION
!!     vte : ELECTRON THERMAL SPEED
!!     fpe : ELECTRON PLASMA FREQUENCY
!!   rmass : ION-TO-ELECTRON MASS RATIO
!!   sigma : (WGE/WPE)**2
!!   beta? : PLASMA BETA FOR ION AND ELECTRON
!!   Ma    : MACH NUMBER OF THE FLOW (IN THE SIMULATION FRAME)
!!   theta : UPSTREAM MAGNETIC FIELD ANGLE NORMAL TO THE X AXIS
!!   lbuf  : INITIAL BUFFER WITHIN WHICH FLOW SMOOTHLY DECREASES TO ZERO
  integer, parameter :: n0    = 64
  real(8), parameter :: vte   = 1.0D0
  real(8), parameter :: fpe   = 1.0D0
  real(8), parameter :: rmass = 25.0D0
  real(8), parameter :: sigma = 0.04D0
  real(8), parameter :: betae = 0.5D0, betai = 0.125D0
  real(8), parameter :: ma    = 3.0
  real(8), parameter :: theta = 90.0/180.0*pi
  real(8), parameter :: lbuf  = 500.0

end module

