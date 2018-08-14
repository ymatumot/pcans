module init

  use const
  use mpi_set
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, allocatable, public :: np2(:,:)
  real(8),              public :: delt
  real(8),              public :: q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
 
  real(8) :: v0, u0, b0, vti, vte, delv


contains

  
  subroutine init__set_param

    use fio, only      : fio__init, fio__input, fio__param
    use boundary, only : boundary__init
    use field, only    : field__init
    use particle, only : particle__init

    real(8) :: fgi, fpi, va, fpe, fge, rgi, rge, ldb
    character(len=64) :: fname_restart

!************** MPI settings  *******************!
    call mpi_set__init(nxgs,nxge,nygs,nyge,nproc)
!*********** End of MPI settings  ***************!
!*********** MEMORY ALLOCATIONS *****************!
    allocate(np2(nys:nye,nsp))
    allocate(uf(6,nxs-2:nxe+2,nys-2:nye+2))
    allocate(up(5,np,nys:nye,nsp))
    allocate(gp(5,np,nys:nye,nsp))
!*********** End of ALLOCATIONS  ***************!

!**** SETTING OTHER NUMERICAL & PHYSICAL CONSTANTS ****!
    r(2) = 1.0D0      ! ELECTRON MASS
    r(1) = r(2)*mr    ! ION MASS
    delt = cfl*delx/c ! TIME STEP SIZE
    ldb  = delx*rdbl
    fpe  = sqrt(beta*rtemp)*c/(sqrt(2.D0)*alpha*ldb)
    fge  = fpe/alpha
    fgi  = fge*r(2)/r(1)
    fpi  = fpe*sqrt(r(2)/r(1))
    va   = fge/fpe*c*sqrt(r(2)/r(1))
    rge  = alpha*ldb*sqrt(2.D0)
    rgi  = rge*sqrt(r(1)/r(2))/sqrt(rtemp)
    vte  = rge*fge
    vti  = vte*sqrt(r(2)/r(1))/sqrt(rtemp)
   
    np2(nys:nye,1) = n0*(nxe+bc-nxs+1)
    np2(nys:nye,2) = np2(nys:nye,1)

    !ELEMENTARY CHARGE
    q(1) = fpi*sqrt(r(1)/(4.0*pi*n0/delx**2))
    q(2) = -q(1)

    !SPEED
    v0 = mf*sqrt(1.D0+5./6.*beta*(1.+rtemp))*va
    u0 = v0/sqrt(1.0D0-(v0/c)**2)

    !HALF SHEAR WIDTH
    delv = dldgi*rgi

    !MAGNETIC FIELD STRENGTH
    b0 = fgi*r(1)*c/q(1)

    !INITIALIZATION OF SUBROUTINES
    call boundary__init(np,nsp,&
                        nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,bc, &
                        nup,ndown,mnpi,mnpr,ncomw,nerr,nstat,   &
                        u0y=(/-0.5*u0, 0.5*u0/))
    call particle__init(np,nsp,&
                        nxs,nxe,nys,nye,nsfo, &
                        q,r,c,delx,delt)
    call field__init(np,nsp,&
                     nxs,nxe,nys,nye,nsfo,bc, &
                     q,c,delx,delt,gfac,      &
                     nup,ndown,mnpr,opsum,ncomw,nerr,nstat)
    call fio__init(np,nsp,&
                   nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsfo,bc, &
                   q,r,c,delx,delt,gfac,n0,                     &
                   fname_param,fname_energy,dir,                &
                   nproc,nrank,nroot,nup,ndown,mnpr,opsum,ncomw,nerr,nstat)
    call random_gen__init(nrank)

    !START FROM RESTART DATA IF AVAILABLE
    if(it0 /= 0)then
       !START FROM THE PAST CALCULATION
       write(fname_restart,'(i6.6,a,i3.3,a)') it0,'_rank=',nrank,'.dat'
       call fio__input(up,uf,np2,it0,fname_restart)
       return
    endif

    !SETUP FOR PARTICLES & FIELDS
    call init__loading
    call fio__param(np2,0.5*r(1)*vti**2,rtemp,fpe,fge,ldb)

  end subroutine init__set_param


  subroutine init__loading

    integer, parameter   :: nbin=nx*20, nchart=nbin, nbin_min=nbin/20
    integer              :: i, j, k, l, ii, isp, isn
    real(8), parameter   :: eps=1d-2
    real(8)              :: dn, dnmax
    real(8)              :: sd, r1, r2, left, right, dx, intf0
    real(8)              :: e0, tp, vfnc, efnc, bfnc, dfnc, gamp, tmp, vd
    real(8)              :: intf(nchart), cx(nchart)
    real(8)              :: bfabs, cdf_n, dcn, vy, gam0, dxvy, dvx, dvy, den, x, y

    !*** FUNCTIONS ***!
    !MAGNETIC FIELD STRENGTH
    bfabs(x) = +0.5*b0*((1.+br)+(1.-br)*tanh((x-x0)/delv))

    !VELOCITY PROFILE
    vy(x) = +0.5*u0*tanh((x-x0)/delv)

    !LORENTZ FACTOR OF MEAN FLOW VELOCITY
    gam0(x) = 1./sqrt(1.D0-(vy(x)/c)**2)

    !GRADIENT OF VELOCITY PROFILE
    dxvy(x) = +0.5*u0/delv*cosh((x-x0)/delv)**(-2)

    !PERTURBATION IN VX
    dvx(x,y) = -0.05*u0*sin(2.*pi*(y-nygs)/(nyge-nygs+1))/cosh((x-x0)/delv)**2

    !NUMBER DENSITY PROFILE
    den(x) = 0.5*n0*((1.+dr)+(1.-dr)*tanh((x-x0)/delv))

    !CDF OF DENSITY PROFILE
    cdf_n(x) = (+((1.+dr)*(   x)+(1.-dr)*delv*log(cosh((   x-x0)/delv))) &
                -((1.+dr)*(nxgs)+(1.-dr)*delv*log(cosh((nxgs-x0)/delv))) &
               )/((1.+dr)*(nxge-nxgs))  

    !CHARGE DENSITY PROFILE
    dcn(x) = abs(bfabs(x)*cos(theta_b)*gam0(x)**2/(4.*pi*q(1)*c)*dxvy(x))

    !drift velocity of ion
    dvy(x) = -c/(4.*pi*q(1)*den(x)*(rtemp+1.)) &
            *0.5*b0*cos(theta_b)/delv*(1.-br)/cosh((x-x0)/delv)**2

    !*** END OF ***!

    !*** SETTING OF FIELDS ***!
    do j=nys-2,nye+2
    do i=nxs-2,nxe+2
       uf(1,i,j) = 0.0D0
       uf(2,i,j) = bfabs(i*delx)*cos(theta_b)
       uf(3,i,j) = bfabs(i*delx)*sin(theta_b)
       uf(4,i,j) = -vy(i*delx)*uf(3,i,j)
       uf(5,i,j) = 0.0
       uf(6,i,j) = 0.0
    enddo
    enddo

    !ADDITIONAL CHARGE TO SATISFY DIV(-VY X BZ) = 4\PI(NI-NE)
    dn = (-uf(4,nxgs,nys)+uf(4,nxge,nys))/(4.*pi*q(1))
    dnmax = 0.0
    do i=nxgs,nxge
       dnmax = max(dnmax,abs(dcn(i*delx)))
    enddo
    dnmax = 0.5*dnmax

    !PARTICLE POSITION - BACKGROUND PARTICLES
    do j=nys,nye
       do ii=1,np2(j,1)-int(dnmax*(nxge+bc-nxgs+1)*delx)
          call random_number(r1)
          left  = nxgs*delx
          right = nxge*delx
          do while(right-left > eps)
             dx = (right+left)/2.
             if(cdf_n(dx) >= r1)then
                right = dx
             else
                left = dx
             endif
          enddo
          !x position
          up(1,ii,j,1) = dx
          up(1,ii,j,2) = up(1,ii,j,1)
          
          !y position
          call random_number(r1)
          up(2,ii,j,1) = dble(j)*delx+r1*delx
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo

    do isp=1,nsp
       isn = sign(1.0D0,dn)
       if(isp == 2) isn = -isn

       !GENERATION OF RANDOM NUMBERS FOR A CHARGE DENSITY PROFILE;
       !CDF IS EVALUATED BY SIMPSON'S RULE
       dx = (nxge+bc-nxgs+1)*delx/nbin
       intf0 = isn*0.5*dcn(nxgs*delx)+dnmax
       do k=1,nbin-2,2
          intf0 = intf0+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                       +2.*(isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax)
       enddo
       k=nbin-1
       intf0 = intf0+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                    +isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax
       intf0 = intf0*dx/3.D0

       do l=1,nchart
          cx(l) = nxgs*delx+(nxge+bc-nxgs+1)*delx*(l-1)/(nchart-1)
          dx = (cx(l)-nxgs*delx)/(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))
          intf(l) = isn*0.5*dcn(nxgs*delx)+dnmax
          do k=1,int(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))-2,2
             intf(l) = intf(l)+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                              +2.*(isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax)
          enddo
          k=(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))-1
          intf(l) = intf(l)+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                           +isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax
          intf(l) = intf(l)*dx/3.D0
       enddo
       intf(1:nchart) = intf(1:nchart)/intf0

       do j=nys,nye
          do ii=np2(j,isp)-int(dnmax*(nxge+bc-nxgs+1)*delx)+1,int(np2(j,isp)-isn*dn/2)
             call random_number(r1)
             left  = 1
             right = nchart
             do while(right-left > 2.)
                if(intf(int((left+right)/2.+0.5)) > r1)then
                   right = (left+right)/2.
                else
                   left = (left+right)/2.
                endif
             enddo
             up(1,ii,j,isp) = cx(int((left+right)/2.+0.5))

             call random_number(r1)
             up(2,ii,j,isp) = dble(j)*delx+delx*r1
          enddo
       enddo
       np2(nys:nye,isp) = np2(nys:nye,isp)-isn*dn/2.
    enddo

    !VELOCITY
    !TOTAL PRESSURE AT X=NXGE
    e0 = -0.5*v0*b0/c
    tp = 0.5*n0*(r(1)*vti**2+r(2)*vte**2)+(b0*b0-e0*e0)/(8.*pi)

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             vfnc = vy(up(1,ii,j,isp))
             dfnc = den(up(1,ii,j,isp))
             bfnc = bfabs(up(1,ii,j,isp))
             efnc = -bfnc*vfnc/c
             
             sd = sqrt( 2.*(tp-(bfnc**2-efnc**2)/(8.*pi)) &
                        /(dfnc*r(1)*(1.+rtemp)) )/sqrt(2.D0)

             if(isp == 2) then
                sd = sd*sqrt(r(1)/r(2)*rtemp)
             endif

             !NON-RELATIVISTIC MAXWELLIAN
             call random_gen__bm(r1,r2)
             up(3,ii,j,isp) = sd*r1
             up(4,ii,j,isp) = sd*r2

             call random_gen__bm(r1,r2)
             up(5,ii,j,isp) = sd*r1

             gamp = sqrt(1.D0+(up(3,ii,j,isp)**2+up(4,ii,j,isp)**2+up(5,ii,j,isp)**2)/(c*c))
             vd = dvy(up(1,ii,j,isp))
             if(isp == 2) vd = -vd*rtemp

             ! DENSITY FIX: ZENITANI, PHYS. PLASMAS 22, 042116 (2015)
             call random_number(r1)
             if(up(3,ii,j,isp)*dvx(up(1,ii,j,isp),up(2,ii,j,isp)) >= 0.)then
                up(3,ii,j,isp) = up(3,ii,j,isp)+dvx(up(1,ii,j,isp),up(2,ii,j,isp))*gamp
             else
                if(r1 < (-dvx(up(1,ii,j,isp),up(2,ii,j,isp))*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = -up(3,ii,j,isp)+dvx(up(1,ii,j,isp),up(2,ii,j,isp))*gamp
                else
                   up(3,ii,j,isp) = up(3,ii,j,isp)+dvx(up(1,ii,j,isp),up(2,ii,j,isp))*gamp
                endif
             endif
             ! DENSITY FIX: ZENITANI, PHYS. PLASMAS 22, 042116 (2015)
             call random_number(r1)
             if(up(4,ii,j,isp)*(+vy(up(1,ii,j,isp))+vd) >= 0.)then
                up(4,ii,j,isp) = (+up(4,ii,j,isp)+(+vy(up(1,ii,j,isp))+vd)*gamp) &
                                *gam0(up(1,ii,j,isp))
             else
                if(r1 < (-(+vy(up(1,ii,j,isp))+vd)*up(4,ii,j,isp)/gamp))then
                   up(4,ii,j,isp) = (-up(4,ii,j,isp)+(+vy(up(1,ii,j,isp))+vd)*gamp) &
                                   *gam0(up(1,ii,j,isp))
                else
                   up(4,ii,j,isp) = (+up(4,ii,j,isp)+(+vy(up(1,ii,j,isp))+vd)*gamp) &
                                   *gam0(up(1,ii,j,isp))
                endif
             endif

          enddo
       enddo
    enddo

  end subroutine init__loading


end module init
