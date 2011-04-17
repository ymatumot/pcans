module init

  use const

  implicit none

  private

  public :: init__set_param, init__inject

  integer, public :: np2(1:nx+bc,nsp)
  integer, public :: itmax, it0, intvl1, intvl2, intvl3
  real(8), public :: delx, delt, gfac
  real(8), public :: c
  real(8), public :: uf(6,0:nx+1)
  real(8), public :: up(4,np,1:nx+bc,nsp)
  real(8), public :: q(nsp), r(nsp), vti, vte, va, rtemp, fpe, fge, rgi, rge, ldb, b0
  !gx, gv, are temporal spaces used for the time integration
  real(8), public :: gp(4,np,1:nx,nsp) !just for initialization
  character(len=6),  public :: dir
  character(len=10), public :: file10
  character(len=10), public :: file12
  integer                   :: n0
  real(8)                   :: pi, v0, u0


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param

    real(8) :: fgi, fpi, alpha, beta
    character(len=14) :: file9 
    character(len=21) :: file11

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
    itmax  = 10000
    intvl1 = 100
    intvl2 = 100
    dir    = './dat/'
    file9  = 'init_param.dat'
    file10 = 'file10.dat'
    file12 = 'energy.dat'
    gfac   = 0.505

    it0    = 0
    if(it0 /= 0)then
       !start from the past calculation
       file11 = '005000_test10.dat'
       call fio__input(up,uf,np2,c,q,r,delt,delx,it0,np,nx,nsp,bc,dir,file11)
       return
    endif

!*********************************************************************
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
!*********************************************************************
    pi   = 4.0*atan(1.0)
    delx = 1.0
    c    = 1.0
    delt = 1.0
    ldb  = delx

    r(1) = 16.0
    r(2) = 1.0

    alpha = 4.0
    beta  = 0.01
    rtemp = 1.0

    fpe = dsqrt(beta*rtemp)*c/(dsqrt(2.D0)*alpha*ldb)
    fge = fpe/alpha

    va  = fge/fpe*c*dsqrt(r(2)/r(1))
    rge = fpe/fge*ldb*dsqrt(2.D0)
    rgi = rge*dsqrt(r(1)/r(2))/dsqrt(rtemp)

    vte = rge*fge
    vti = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)

    fgi = fge*r(2)/r(1)
    fpi = fpe*dsqrt(r(2)/r(1))

    n0  = 100
    np2(1:nx+bc,1:2) = n0
    
    if(max(np2(1,1), np2(nx+bc,1), np) > np)then
       write(*,*)'Too large number of particles'
       stop
    endif

    !charge
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*n0))
    q(2) = -q(1)

    !Magnetic field strength
    b0 = fgi*r(1)*c/q(1)

    call init__loading
    call init__set_field
    call fio__param(np,nx,nsp,np2,c,q,r,vti,vte,va,rtemp,fpe,fge,ldb,delt,delx,bc,dir,file9)

  end subroutine init__set_param


  subroutine init__loading

    use boundary, only : boundary__particle

    integer :: i, ii, isp
    real(8) :: sd, aa, bb

    call random_seed()

    !particle position
    isp = 1
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          call random_number(aa)
          up(1,ii,i,1) = dble(i)+delx*aa
          up(1,ii,i,2) = up(1,ii,i,1) 
       enddo
    enddo

    v0 = 0.1*c
    u0 = v0/dsqrt(1.-(v0/c)**2)

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/dsqrt(2.0D0)
       endif
       if(isp == 2) then
          sd = vte/dsqrt(2.0D0)
       endif

       do i=1,nx+bc
          do ii=1,np2(i,isp)
             call random_number(aa)
             call random_number(bb)
             up(2,ii,i,isp) = sd*dsqrt(-2.*dlog(aa))*cos(2.*pi*bb)+u0
             up(3,ii,i,isp) = sd*dsqrt(-2.*dlog(aa))*sin(2.*pi*bb)

             call random_number(aa)
             call random_number(bb)
             up(4,ii,i,isp) = sd*dsqrt(-2.*dlog(aa))*cos(2.*pi*bb)
          enddo
       enddo
    enddo

    call boundary__particle(up,np,nx,nsp,np2,bc)

  end subroutine init__loading


  subroutine init__set_field

    use boundary, only : boundary__field

    integer :: i

    !magnetic field
    do i=1,nx+bc
       uf(1,i) = 0.0
    enddo
    do i=1,nx
       uf(2,i) = 0.0
       uf(3,i) = b0
    enddo

    !electric field
    do i=1,nx
       uf(4,i) = 0.0
    enddo
    do i=1,nx+bc
       uf(5,i) = v0*b0/c
       uf(6,i) = 0.0
    enddo

    call boundary__field(uf,nx,bc)

  end subroutine init__set_field


  subroutine init__inject

    use boundary, only : boundary__particle

    integer :: isp, i, ii, ii2, ii3, dn
    real(8) :: sd, aa, bb
    !Inject particles at x=1

    i   = 1
    isp = 1
    dn  = n0-max(np2(i,1),np2(i,2))

    do ii=1,dn
       ii2 = np2(i,1)+ii
       ii3 = np2(i,2)+ii
       call random_number(aa)
       up(1,ii2,i,1) = dble(i)+delx*aa
       up(1,ii3,i,2) = up(1,ii2,i,1)
    enddo

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp == 1) then 
          sd = vti/dsqrt(2.0D0)
       endif
       if(isp == 2) then
          sd = vte/dsqrt(2.0D0)
       endif

       do ii=np2(i,isp)+1,n0
          call random_number(aa)
          call random_number(bb)
          up(2,ii,i,isp) = sd*dsqrt(-2.*dlog(aa))*cos(2.*pi*bb)+u0
          up(3,ii,i,isp) = sd*dsqrt(-2.*dlog(aa))*sin(2.*pi*bb)

          call random_number(aa)
          call random_number(bb)
          up(4,ii,i,isp) = sd*dsqrt(-2.*dlog(aa))*cos(2.*pi*bb)
       enddo
    enddo

    do isp=1,nsp
       np2(i,isp) = np2(i,isp)+max(dn,0)
    enddo

    !set Ex and Bz
    i=1
    uf(3,i) = b0
    uf(5,i) = v0*b0/c

    call boundary__particle(up,np,nx,nsp,np2,bc)

  end subroutine init__inject

end module init
