module init

  use const
  use mpi_set

  implicit none

  private

  public :: init__set_param

  integer, public, parameter :: nroot=0
  integer, allocatable, public :: np2(:,:)
  integer, public :: itmax, it0, intvl1, intvl2, intvl3
  real(8), public :: delx, delt, gfac
  real(8), public :: c
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), public :: q(nsp), r(nsp)
  real(8), allocatable, public :: gp(:,:,:,:)
  character(len=64), public :: dir
  character(len=64), public :: file12
  real(8)                   :: pi, vti, vte, va, rtemp, fpe, fge, rgi, rge, ldb, b0
  real(8)                   :: v0, Ti, Te, nrate, lcs
  integer                   :: ncs, nbg


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param
    real(8) :: fgi, fpi, alpha, beta, n0
    character(len=64) :: file9 
    character(len=64) :: file11
    integer :: ii, j
!************** MPI settings  *******************!
    call mpi_set__init(nxgs,nxge,nygs,nyge,nproc)

    allocate(np2(nys:nye,nsp))
    allocate(uf(6,nxs1:nxe1,nys1:nye1))
    allocate(up(5,np,nys:nye,nsp))
    allocate(gp(5,np,nys:nye,nsp))

!*********** End of MPI settings  ***************!

!*********************************************************************
!   time0   : start time (if time0 < 0, initial data from input.f)
!   itmax   : number of iteration
!   it0     : base count
!   intvl1  : storage interval for particles & fields
!   intvl2  : printing interval for energy variation
!   intvl3  : printing interval for wave analysis
!   dir     : directory name for data output
!   file??  : output file name for unit number ??
!           :  9 - initial parameters
!           : 10 - for saving all data
!           : 11 - for starting from saved data
!           : 12 - for saving energy history
!   gfac    : implicit factor
!             gfac < 0.5 : unstable
!             gfac = 0.5 : no implicit
!             gfac = 1.0 : full implicit
!*********************************************************************
    pi     = 4.0*atan(1.0)
    itmax  = 3000

    intvl1 = 1000 ! full data
    intvl2 = 100 ! energy
    intvl3 = 1000000 ! wave analysis

    dir    = './dat/'
    file9  = 'init_param.dat'
    file12 = 'energy.dat'
    gfac   = 0.505
    it0    = 0



    it0    = 0
    if(it0 /= 0)then
       !start from the past calculation
       write(file11,'(a,i3.3,a)')'000100_rank=',nrank,'.dat'
       call fio__input(up,uf,np2,c,q,r,delt,delx,it0,                             &
                       np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,bc,nproc,nrank, &
                       dir,file11)
       return
    endif


!*********************************************************************
!  NOTE: 
!    following parameters are defined in the framework of the non-
!   relativistic plasma. The relativistic collection is needed
!   to see the actual value (e.g., wce -> wce / \gamma). 
!
!   
!   r(1)  : ion mass             r(2)  : electron mass
!   q(1)  : ion charge           q(2)  : electron charge
!   c     : speed of light       ldb   : debye length
!
!   rgi   : ion Larmor radius    rge   : electron Larmor radius
!   fgi   : ion gyro-frequency   fge   : electron gyro-frequency
!   vti   : ion thermal speed    vte   : electron thermal speed
!   b0    : magnetic field       
!  
!   alpha : wpe/wge
!   beta  : ion plasma beta
!   rtemp : Te/Ti
!   lcs   : thickness of the current sheet
!   n0    : number of super-particles in the cell inside the current sheet
!   nrate : ratio of the particle number density between the 
!           background plasma and the current sheet
!*********************************************************************
    delx  = 1.0
    c     = 1.0
    delt  = 0.2*delx/c

    ldb   = delx*3.d0
    r(1)  = 1.0
    r(2)  = 1.0

    alpha = 0.5d0
    beta  = 0.5d0
    rtemp = 1.d0

    n0    = 100
    nrate = 0.1d0
    lcs   = delx*10



    ! drift velocity
    v0 = ldb/lcs

    ! number of super-particles
    ncs = n0
    nbg = n0*nrate

    fpe = dsqrt(beta*rtemp)*c/(dsqrt(2.D0)*alpha*ldb)
    fge = fpe/alpha

    va  = fge/fpe*c*dsqrt(r(2)/r(1))
    rge = fpe/fge*ldb*dsqrt(2.D0)
    rgi = rge*dsqrt(r(1)/r(2))/dsqrt(rtemp)

    vte = rge*fge
    vti = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)

    fgi = fge*r(2)/r(1)
    fpi = fpe*dsqrt(r(2)/r(1))

    if(nrank == nroot)then
       if(max(np2(nys,1), np2(nye,1), np) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif

    ! charge
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*n0))
    q(2) = -q(1)

    ! Magnetic field strength
    b0 = fgi*r(1)*c/q(1)

    np2(nys:nye,1:nsp) = nbg*(nxe+bc-nxs+1) + 2*lcs*ncs*2


    ! temperature
    Ti = beta*B0*B0/(8.d0*pi*ncs)
    Te = rtemp*Ti

    call init__loading
    call fio__param(np,nsp,np2,                     &
                    nxgs,nxge,nygs,nyge,nys,nye,    &
                    c,q,r,n0,0.5*r(1)*vti**2,rtemp,fpe,fge, &
                    ldb,delt,delx,dir,file9,        &
                    nroot,nrank)
  end subroutine init__set_param



  subroutine init__loading

    use boundary, only : boundary__particle, boundary__field
    integer, parameter   :: nbin=nx*20, nchart=nbin, nbin_min=nbin/20
    integer              :: i, j, k, l, ii, isp
    real(8), parameter   :: eps=1d-2
    real(8)              :: aa, x, y, xs, ys
    real(8) :: ih, jh

    ! position of the reconnection point
    real(8) :: x1, y1, x2, y2, x3, y3, xh, yh, xu, xd
    real(8) :: rr1, rr2, rr3,delb
    integer :: ixu, ixd, ixh, iyh

    ! random number
    integer,dimension(nys:nye) :: iis,iie
    integer :: iran, nptot, ir
    integer :: nran
    real(8),dimension(:),allocatable :: ran
    real(8) :: norm,Tem,vbulk,umax,uc,fc,uh,fh,w,wl,wr,ifc,ifh,gamb,bb,fstep,fu,csth,snth,wkuz,wkut,uwk,ewk,pwk,phi,ua


    ! vector potential defined at cell center
    real(8),dimension(:,:),allocatable :: az
    
    allocate(az(nxs1-1:nxe1,nys1-1:nye1))

    ! set initial parameter
    ! reconnection point is set at
    ! (x1, y1), (x2, y2), and (x3, y3)
    ixh = (nxge - nxgs)*0.5d0  ; xh = ixh*delx
    iyh = (nyge - nygs)*0.5d0  ; yh = iyh*delx
    ixu = (nxge - nxgs)*0.75d0 ; xu = ixu*delx
    ixd = (nxge - nxgs)*0.25d0 ; xd = ixd*delx
    x1  = xd; y1 = DBLE(nygs)*delx
    x2  = xd; y2 = DBLE(nyge)*delx
    x3  = xu; y3 = yh
    
    !*** setting of fields ***!
    !* Magnetic Field
    ! magnetic field is obtained by descretizing the vector potential.
    ! delb is the relative amplitude of the perturbed field.
    delb=0.1d0
    do j = nys1-1, nye1
    do i = nxs1-1, nxe1
       xs  = dble(i+0.5)*delx
       ys  = dble(j+0.5)*delx

       rr1 = ((xs - x1)**2.d0 + (ys - y1)**2.d0)/(lcs**2.d0)
       rr2 = ((xs - x2)**2.d0 + (ys - y2)**2.d0)/(lcs**2.d0)
       rr3 = ((xs - x3)**2.d0 + (ys - y3)**2.d0)/(lcs**2.d0)
       az(i,j) =-b0*(xs &
                    - lcs*log(cosh((xs-xd)/lcs))  &
                    + lcs*log(cosh((xs-xu)/lcs)) &
                    - 2.d0*lcs*delb*exp(-rr1*0.25d0) &
                    - 2.d0*lcs*delb*exp(-rr2*0.25d0) &
                    + 2.d0*lcs*delb*exp(-rr3*0.25d0))
    enddo
    enddo
    do j=nys1,nye1
    do i=nxs1,nxe1
       uf(1,i,j) = (az(i,j) - az(i,j-1))/delx
       uf(2,i,j) =-(az(i,j) - az(i-1,j))/delx
       uf(3,i,j) = 0.d0
    enddo
    enddo
    deallocate(az)

    !* Electric Field
    do j=nys,nye
    do i=nxs,nxe
       uf(4,i,j) = 0.d0
       uf(5,i,j) = 0.d0
       uf(6,i,j) = 0.d0
    enddo
    enddo


    call boundary__field(uf,                 &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)
    !*** end of ***!
!---------------------------------------------------------------

    !*** setting of particles ***!
    ! Plasma is consist of background plasma and the
    ! harris sheet. 
    ! There are two current sheets in the simulation box
    ! because the double periodic boundary condition is adopted.
    ! Particle distribution function is given by the
    ! relativistic Maxwellian.

    ! generate random number
    iran = 1
    call random__init(iran,nrank)
    do j = nys, nye
    do ii =1, np2(j,1)
       nptot = nptot + 1
    enddo
    enddo
    nran = nptot*2
    allocate(ran(nran))
    call random__generate(nran,ran,iran)    


    ! *** particle position 
    ir =1
    isp=1
    do j=nys,nye
       ! background particles
       do ii=1, nbg*(nxe+bc-nxs+1) 
          !x position
          up(1,ii,j,1) = nxs*delx + ran(ir)*delx*(nxe+bc-nxs+1.d0)
          up(1,ii,j,2) = up(1,ii,j,1)
          
          !y position
          up(2,ii,j,1) = dble(j)*delx + ran(ir+1)*delx
          up(2,ii,j,2) = up(2,ii,j,1)
          
          ir = ir + 2
       enddo

       ! 1st Harris sheet
       do ii = nbg*(nxe+bc-nxs+1)+1, nbg*(nxe+bc-nxs+1)+2*lcs*ncs
          aa=2.d0*ran(ir)-1.d0
          up(1,ii,j,1) = lcs*log((1 + aa)/(1.d0-aa))*0.5d0 + xd
          up(1,ii,j,2) = up(1,ii,j,1)

          up(2,ii,j,1) = dble(j)*delx + ran(ir+1)*delx
          up(2,ii,j,2) = up(2,ii,j,1)

          ir = ir + 2
       enddo
       ! 2nd Harris sheet
       do ii = nbg*(nxe+bc-nxs+1)+2*lcs*ncs+1, np2(j,1)
          aa=2.d0*ran(ir)-1.d0
          up(1,ii,j,1) = lcs*log((1 + aa)/(1.d0-aa))*0.5d0 + xu
          up(1,ii,j,2) = up(1,ii,j,1)

          up(2,ii,j,1) = dble(j)*delx + ran(ir+1)*delx
          up(2,ii,j,2) = up(2,ii,j,1)

          ir = ir + 2
       enddo
    enddo
    deallocate(ran)

    ! *** particle velocity
    ! background plasma
    iis = 1; iie = int(nbg*(nxe+bc-nxs-1))
    call thermal__r(up,5,np,nys,nye,nsp,1,iis,iie,nys,nye,Ti,r(1),&
         Ti*100.d0,c,iran)
    call thermal__r(up,5,np,nys,nye,nsp,2,iis,iie,nys,nye,Te,r(2),&
         Te*100.d0,c,iran)

    ! 1st current sheet
    iis = iie+1; iie = int(nbg*(nxe+bc-nxs+1))+2*lcs*ncs
    call thermal__r_shift(up,5,np,nys,nye,nsp,1,iis,iie,nys,nye,&
         Ti,r(1),-v0,3,Ti*100.d0,c,iran)
    call thermal__r_shift(up,5,np,nys,nye,nsp,2,iis,iie,nys,nye,&
         Te,r(2), v0,3,Te*100.d0,c,iran)

    ! 2nd current sheet
    iis = iie+1;  iie(:) = np2(:,1)
    call thermal__r_shift(up,5,np,nys,nye,nsp,1,iis,iie,nys,nye,&
         Ti,r(1), v0,3,Ti*100.d0,c,iran)
    call thermal__r_shift(up,5,np,nys,nye,nsp,2,iis,iie,nys,nye,&
         Te,r(2),-v0,3,Te*100.d0,c,iran)


    call boundary__particle(up,                                        &
                            np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye,bc, &
                            nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

  end subroutine init__loading


end module init
