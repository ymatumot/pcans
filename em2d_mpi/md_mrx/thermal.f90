!========================================================|
! thermal.f90
!
!  PURPOSE: make maxwell distribution function 
!  This file is consist of ... 
!
!   thermal_r:
!    relativistic maxwellian is generated using reduction
!    method
!
!   thermal_r_shift:
!    relativistic shifted maxwellian is generated using 
!    reduction method
! 
!                                           H. Takahashi
!========================================================|
subroutine thermal__r(up,nv,np,nys,nye,nsp,isp,iis,iie,js,je,T,m,id,umax,c,iran)
  implicit none

  integer,intent(in)               :: nv, np,nys, nye, nsp
  integer,intent(in)               :: isp, js, je,iran
  integer,intent(in),dimension(nys:nye) :: iis, iie
  integer,intent(in)               :: id
  real(8),intent(in)               :: T, m, umax, c
  real(8),intent(inout),dimension(nv,np,nys:nye,nsp) :: up

  real(8) :: Tem
  real(8) :: uc, fc, uh, fh, wl, wr, w, fu, fstep
  real(8) :: ua, csth, snth, phi, norm
  real(8) :: aa, bb
  real(8) :: ifc, ifh
  real(8) :: pi  
  integer :: ivs
  integer :: j, ii
  real(8),allocatable :: ran(:)
  integer :: nran, ir


  ivs = nv-2
  pi  = 4.d0*atan(1.d0)
  
  nran = 0
  do j = js, je
  do ii = iis(j), iie(j)
     nran = nran + 1
  enddo
  enddo

  nran = 2*nran*21
  allocate(ran(nran))
  call random__generate(nran,ran,iran)
  norm = 1.0D+3
  Tem = T/m

  ! peak value of the distribution function
  uc =sqrt(2.d0*Tem*Tem*(1.d0 + sqrt(1.d0+1.d0/(Tem*Tem))))
  fc =norm*exp(-sqrt(uc*uc+1.d0)/Tem)*uc*uc

  ! make a step function
  uh =(uc + umax)*0.5d0
  fh =norm*exp(-sqrt(uh*uh+1.d0)/Tem)*uh*uh


  ! weigth
  wl = fc*uh
  wr = fh*(umax-uh)
  w  = 1.d0/(wl+wr)
  wl = wl*w
  wr = wr*w
  ifc= uh/wl
  ifh= (umax-uh)/(1.d0-wl)

  ir = 1
  do j = js, je
  do ii = iis(j), iie(j)
     fu = 1.d0; fstep = 10.d0
     do while(fu .lt. fstep)
        if(ran(ir).lt. wl)then
           ! left hand of the step function
           ua = ran(ir)*ifc
           fstep = fc*ran(ir+1)
        else
           ! right hand of the step function
           ua = (ran(ir) - wl)*ifh + uh
           fstep = fh*ran(ir+1)
        end if
        ir = ir + 2
        fu = norm*exp(-sqrt(ua*ua + 1.d0)/Tem)*ua*ua
        if(ir .ge. nran-2)then
           call random__generate(nran,ran,iran)
           ir = 1
        endif

     end do

     csth = 1.d0-2.d0*ran(ir)
     snth = sqrt(1.d0 - csth*csth)
     phi  = 2.d0*pi*ran(ir+1)
     uc   = ua*c
     ir   = ir + 2

     up(ivs,ii,j,isp)   = uc*snth*cos(phi)
     up(ivs+1,ii,j,isp) = uc*snth*sin(phi)
     up(ivs+2,ii,j,isp) = uc*csth
  enddo
  enddo

  deallocate(ran)
  return
end subroutine thermal__r



subroutine thermal__r_shift(up,nv,np,nys,nye,nsp,isp,iis,iie,js,je,T,m,v0,idir,id,umax,c,iran)
  implicit none

  integer,intent(in)               :: nv, np,nys, nye, nsp
  integer,intent(in)               :: isp, js, je, idir, iran
  integer,intent(in),dimension(nys:nye) :: iis, iie
  integer,intent(in)               :: id
  real(8),intent(in)               :: T, m, v0, umax, c
  real(8),intent(inout),dimension(nv,np,nys:nye,nsp) :: up

  real(8) :: Tem
  real(8) :: uc, fc, uh, fh, wl, wr, w, fu, fstep
  real(8) :: ua, csth, snth, phi, norm
  real(8) :: aa, bb
  real(8) :: ifc, ifh
  real(8) :: gamb, vb,wkuz,wkut, uwk, pwk, ewk
  real(8) :: pi  
  integer :: ivs
  integer :: j, ii
  integer :: id1, id2, id3
  integer :: nseed
  integer,allocatable :: seed(:)
  real(8),allocatable :: ran(:)
  integer :: nran, ir
  integer :: ic
  integer :: tot, mis
  ivs = nv-2
  pi  = 4.d0*atan(1.d0)
  id3 = ivs
  id1 = ivs + 1
  id2 = ivs + 2

  select case(idir)
  case(1)
     id3 = ivs
     id1 = ivs + 1
     id2 = ivs + 2
  case(2)
     id3 = ivs+1
     id1 = ivs 
     id2 = ivs + 2
  case(3)
     id3 = ivs+2
     id1 = ivs 
     id2 = ivs + 1     
  end select

  nran = 0
  do j = js, je
  do ii = iis(j), iie(j)
     nran = nran + 1
  enddo
  enddo
  nran = 1*nran + 3*nran*10
  allocate(ran(nran))
  call random__generate(nran,ran,iran)

  norm = 1.d+3
  Tem = T/m





  ! peak value of the distribution function
  uc =sqrt(2.d0*Tem*Tem*(1.d0 + sqrt(1.d0 + 1.d0/(Tem*Tem))))
  fc =norm*exp(-sqrt(uc*uc + 1.d0)/Tem)*uc*uc


  ! set step function
  uh  =(uc + umax)*0.5d0
  fh  = norm*exp(-sqrt(uh*uh + 1.d0)/Tem)*uh*uh

  ! weigth
  wl  = fc*uh
  wr  = fh*(umax - uh)

  w   = 1.d0/(wl + wr)
  wl  = wl*w
  wr  = wr*w
  ifc = uh/wl
  ifh = (umax - uh)/(1.d0 - wl)
  vb  = v0/c
  gamb=1.d0/sqrt(1.d0 - vb*vb)

  ir = 1

  tot=0
  mis=0
  do j = js,je
  do ii = iis(j), iie(j)

     fu = 1.d0; fstep = 10.d0
     do while(fu .lt. fstep)
        tot=tot+1
        mis=mis+1
        if(ran(ir).lt. wl)then
           ! left hand of the step function
           ua = ran(ir)*ifc
           fstep = fc*ran(ir+1)
        else
           ! right hand of the step function
           ua = (ran(ir) - wl)*ifh+uh
           fstep = fh*ran(ir+1)
        end if
        fu = norm*exp(-sqrt(ua*ua + 1.d0)/Tem)*ua*ua
        
        csth = 1.d0 - 2.d0*ran(ir+2)
        snth = sqrt(1.d0 - csth*csth)
        
        wkuz = ua*csth
        wkut = sqrt(1.d0 + ua*ua)

        ! Lorentz transformation
        uwk  = (   wkut + vb*wkuz)
        pwk  = (vb*wkut +    wkuz)
        ewk  = uwk/wkut/(gamb*(1.d0 + abs(vb)))
        fu   = fu*ewk

        ir = ir + 3
        if(ir .ge. nran-3)then
           call random__generate(nran,ran,iran)
           ir = 1
        endif
     end do
     mis=mis-1
     phi   = 2.d0*pi*ran(ir)
     uc    = ua*c
     up(id1,ii,j,isp) = uc*snth*cos(phi)
     up(id2,ii,j,isp) = uc*snth*sin(phi)
     up(id3,ii,j,isp) = gamb*(vb*wkut+  wkuz)*c
     ir = ir + 1
  end do
  end do

  return
end subroutine thermal__r_shift
