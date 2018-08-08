module const

  implicit none

!!************************ NUMERICAL CONSTANTS ***********************************!!
  integer, parameter :: nx  = 512   ! NUMBER OF GRID POINTS
  integer, parameter :: np  = 5000  ! NUMBER OF PARTICLES IN EACH CELL
  integer, parameter :: nsp = 2     ! NUMBER OF PARTICLE SPECIES
  integer, parameter :: bc  = 0     ! BOUNDARY CONDITION (PERIODIC:0, REFLECTIVE:-1)

!! SETUP FOR SUBROUTINES CALLED IN MAIN PROGRAM
  integer, parameter :: itmax  = 80000 !NUMBER OF ITERATION
  integer            :: it0    = 0    !0:INITIAL, NONZERO/9999999: RESTART DATA
  integer, parameter :: intvl1 = 500  !INTERVAL FOR PARTICLES & FIELDS STORAGE
  integer, parameter :: intvl2 = 500  !INTERVAL FOR ENERGY CALC.
  integer, parameter :: intvl3 = 100  !INTERVAL FOR RECORDING MOMENT & FIELDS DATA	
  character(len=128) :: dir    = './dat/' !DIRECTORY FOR OUTPUT
  character(len=128) :: dir_mom= './mom/' !DIRECTORY FOR MOMENT DATA
  character(len=128) :: dir_psd= './psd/' !DIRECTORY FOR MOMENT DATA
  character(len=128) :: file9  = 'init_param.dat' !FILE NAME OF INIT CONDITIONS
  character(len=128) :: file10 = 'file10.dat' !FILE NAME OF PARTICLE DATA
  character(len=128) :: file12 = 'energy.dat' !FILE NAME OF ENERGY VARIATION

!! OTHER CONSTANTS
  real(8), parameter :: cfl    = 0.5D0   !CFL CONDITION FOR LIGHT WAVE
  real(8), parameter :: gfac   = 0.501D0 !IMPLICITNESS FACTOR > 0.5
  real(8), parameter :: rdbl   = 1.0D0   !DEBYE LENGTH / CELL WIDTH
  real(8), parameter :: pi     = 4.0D0*atan(1.0D0)

!!************************ PHYSICAL CONSTANTS ***********************************!!
!!      n0 : NUMBER OF PARTICLES/CELL
!!     vte : ELECTRON THERMAL SPEED
!!     fpe : ELECTRON PLASMA FREQUENCY
!!   rmass : ION-TO-ELECTRON MASS RATIO
!!   sigma : (WGE/WPE)**2
!!   beta? : PLASMA BETA FOR ION AND ELECTRON
!!   nbeam : DENSITY RATIO (BEAM / TOTAL)
!!   vbeam : RELATIVE STREAMING VELOCITY
!!   tbeam : TEMPERATURE RATIO (BEAM / CORE)
  integer, parameter :: n0     = 10
  real(8), parameter :: vte    = 1.0D0
  real(8), parameter :: fpe    = 1.0D0
  real(8), parameter :: rmass  = 100.0D0
  real(8), parameter :: sigma  = 0.01D0
  real(8), parameter :: betai  = 0.5D0, betae = 0.005D0
  real(8), parameter :: nbeam  = 0.5, vbeam = 20.0, tbeam = 1.0

end module
