  subroutine field__cip(uf,up,gp,np,nx,nsp,np2,q,c,delx,delt,bc)
!*********************************************
!
! EM field solver by CIP method
!
!*********************************************

    use boundary, only : boundary__field, boundary__curre,  boundary__particle

    integer, intent(in)    :: np, nx, nsp, bc
    integer, intent(in)    :: np2(1:nx+bc,nsp)
    real(8), intent(in)    :: q(nsp), c, delx, delt
    real(8), intent(in)    :: gp(4,np,1:nx+bc,nsp)
    real(8), intent(inout) :: up(4,np,1:nx+bc,nsp)
    real(8), intent(inout) :: uf(6,0:nx+1)
    integer                    :: ii, i, isp, iup
    integer, save              :: iflag
    real(8), parameter         :: eps=1d-3
    real(8), allocatable, save :: duf(:,:)
    real(8)                    :: uj(3,-1:nx+2)
    real(8)                    :: gf(6,0:nx+1)
    real(8)                    :: dgf(6,0:nx+1)
    real(8)                    :: d2, d, di, dl, tmp, tmp2, ab, sx, c1, c2, c3, alpha, beta, dS, fac
    real(8)                    :: pi

    pi = 4.0*atan(1.0)

    !position at n+1/2
    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             up(1,ii,i,isp) = 0.5*(up(1,ii,i,isp)+gp(1,ii,i,isp))
             up(2:4,ii,i,isp) = gp(2:4,ii,i,isp)
          enddo
       enddo
    enddo
    call ele_cur_cip(uj,up,np,nx,nsp,np2,bc,q,c)
    call boundary__curre(uj,nx,bc)

    !convert from E and B to Ex+-By and Ey+-Bx
    uf(5,0:nx+1) = uf(5,0:nx+1)+uf(3,0:nx+1)     !Ey+Bz
    uf(3,0:nx+1) = uf(5,0:nx+1)-2.0*uf(3,0:nx+1) !Ey-Bz
    uf(6,0:nx+1) = uf(6,0:nx+1)+uf(2,0:nx+1)     !Ez+By
    uf(2,0:nx+1) = uf(6,0:nx+1)-2.0*uf(2,0:nx+1) !Ez-By

    !gradient at t=0
    if(iflag /= 1)then
       allocate(duf(6,0:nx+1))
       do i=1,nx
          duf(2,i) = 0.5*(-uf(2,i-1)+uf(2,i+1))
          duf(3,i) = 0.5*(-uf(3,i-1)+uf(3,i+1))
          duf(5,i) = 0.5*(-uf(5,i-1)+uf(5,i+1))
          duf(6,i) = 0.5*(-uf(6,i-1)+uf(6,i+1))
       enddo
       call boundary__field(duf,nx,bc)
       iflag = 1
    endif

    !+non advective term - first step
    fac = +4.*pi*delt*0.5
    uf(4,1:nx) = uf(4,1:nx)-fac*uj(1,1:nx)
    uf(5,1:nx) = uf(5,1:nx)-fac*uj(2,1:nx)
    uf(3,1:nx) = uf(3,1:nx)-fac*uj(2,1:nx)
    uf(6,1:nx) = uf(6,1:nx)-fac*uj(3,1:nx)
    uf(2,1:nx) = uf(2,1:nx)-fac*uj(3,1:nx)
    !update gradient
    call boundary__field(duf,nx,bc)
    do i=1,nx
       dS = -fac*0.5*(-uj(2,i-1)+uj(2,i+1))
       duf(3,i) = duf(3,i)+dS
       duf(5,i) = duf(5,i)+dS
    enddo
    do i=1,nx
       dS = -fac*0.5*(-uj(3,i-1)+uj(3,i+1))
       duf(2,i) = duf(2,i)+dS
       duf(6,i) = duf(6,i)+dS
    enddo

    call boundary__field(uf,nx,bc)
    call boundary__field(duf,nx,bc)
       
    !set constant
    d2  = 1.D0/(delx*delx)

    !Advection of Ey+Bz
    iup = -1
    d  = dble(iup)*delx
    di = 1./d
    dl = -c*delt
    do i=1,nx
       sx = (uf(5,i+iup)-uf(5,i))*di

       if(duf(5,i+iup)*duf(5,i) < 0)then
          alpha = 1.0
       else
          alpha = 0.0
       endif
       
       tmp = duf(5,i+iup)-sx + 1D-10

       beta = 0.5*(1.+tanh((dabs(tmp)-eps)/(0.1*eps))) &
             *(dabs((sx-duf(5,i))/tmp)-1.0)*di
       ab = alpha*beta

       c3 = ( duf(5,i)-sx+tmp*(1.+ab*d) )*d2
       c2 = sx*ab+(sx-duf(5,i))*di-c3*d
       c1 = duf(5,i)+uf(5,i)*ab

       tmp2 = 1./(1.+ab*dl)

       gf(5,i) = (((c3*dl+c2)*dl+c1)*dl+uf(5,i))*tmp2
       dgf(5,i) = ( ((3.*c3*dl+2.*c2)*dl+c1) &
                     -(((c3*dl+c2)*dl+c1)*dl+uf(5,i))*tmp2*ab )*tmp2
    enddo

    !Advection of Ey-Bz
    iup = +1
    d  = dble(iup)*delx
    di = 1./d
    dl = +c*delt
    do i=1,nx
       sx = (uf(3,i+iup)-uf(3,i))*di

       if(duf(3,i+iup)*duf(3,i) < 0)then
          alpha = 1.0
       else
          alpha = 0.0
       endif
       
       tmp = duf(3,i+iup)-sx + 1D-10
       beta = 0.5*(1.+tanh((dabs(tmp)-eps)/(0.1*eps))) &
             *(dabs((sx-duf(3,i))/tmp)-1.0)*di

       ab = alpha*beta
       
       c3 = ( duf(3,i)-sx+tmp*(1.+ab*d) )*d2
       c2 = sx*ab+(sx-duf(3,i))*di-c3*d
       c1 = duf(3,i)+uf(3,i)*ab

       tmp2 = 1./(1.+ab*dl)

        gf(3,i) = (((c3*dl+c2)*dl+c1)*dl+uf(3,i))*tmp2
       dgf(3,i) = ( ((3.*c3*dl+2.*c2)*dl+c1) &
                     -(((c3*dl+c2)*dl+c1)*dl+uf(3,i))*tmp2*ab )*tmp2
    enddo
    !Advection of Ez+By
    iup = +1
    d  = dble(iup)*delx
    di = 1./d
    dl = +c*delt
    do i=1,nx
       sx = (uf(6,i+iup)-uf(6,i))*di

       if(duf(6,i+iup)*duf(6,i) < 0)then
          alpha = 1.0
       else
          alpha = 0.0
       endif
       
       tmp = duf(6,i+iup)-sx + 1D-10
       beta = 0.5*(1.+tanh((dabs(tmp)-eps)/(0.1*eps))) &
             *(dabs((sx-duf(6,i))/tmp)-1.0)*di

       ab = alpha*beta
       
       c3 = ( duf(6,i)-sx+tmp*(1.+ab*d) )*d2
       c2 = sx*ab+(sx-duf(6,i))*di-c3*d
       c1 = duf(6,i)+uf(6,i)*ab

       tmp2 = 1./(1.+ab*dl)

        gf(6,i) = (((c3*dl+c2)*dl+c1)*dl+uf(6,i))*tmp2
       dgf(6,i) = ( ((3.*c3*dl+2.*c2)*dl+c1) &
                     -(((c3*dl+c2)*dl+c1)*dl+uf(6,i))*tmp2*ab )*tmp2
    enddo
    !Advection of Ez-By
    iup = -1
    d  = dble(iup)*delx
    di = 1./d
    dl = -c*delt
    do i=1,nx
       sx = (uf(2,i+iup)-uf(2,i))*di

       if(duf(2,i+iup)*duf(2,i) < 0)then
          alpha = 1.0
       else
          alpha = 0.0
       endif
       
       tmp = duf(2,i+iup)-sx + 1D-10
       beta = 0.5*(1.+tanh((dabs(tmp)-eps)/(0.1*eps))) &
             *(dabs((sx-duf(2,i))/tmp)-1.0)*di

       ab = alpha*beta
       
       c3 = ( duf(2,i)-sx+tmp*(1.+ab*d) )*d2
       c2 = sx*ab+(sx-duf(2,i))*di-c3*d
       c1 = duf(2,i)+uf(2,i)*ab

       tmp2 = 1./(1.+ab*dl)

       gf(2,i) = (((c3*dl+c2)*dl+c1)*dl+uf(2,i))*tmp2
       dgf(2,i) = ( ((3.*c3*dl+2.*c2)*dl+c1) &
                     -(((c3*dl+c2)*dl+c1)*dl+uf(2,i))*tmp2*ab )*tmp2
    enddo

    !update and convert from Ey+-Bz and Ez+-By to E and B
    !+non advective term
    uf(2,1:nx) = 0.5*(-gf(2,1:nx)+gf(6,1:nx))
    uf(3,1:nx) = 0.5*(-gf(3,1:nx)+gf(5,1:nx))
    uf(4,1:nx) = uf(4,1:nx)-fac*uj(1,1:nx)
    uf(5,1:nx) = 0.5*(+gf(3,1:nx)+gf(5,1:nx))-fac*uj(2,1:nx)
    uf(6,1:nx) = 0.5*(+gf(2,1:nx)+gf(6,1:nx))-fac*uj(3,1:nx)

    !update gradient
    !non-advective term
    do i=1,nx
       dS = -fac*0.5*(-uj(2,i-1)+uj(2,i+1))
       duf(3,i) = dgf(3,i)+dS
       duf(5,i) = dgf(5,i)+dS
    enddo

    do i=1,nx
       dS = -fac*0.5*(-uj(3,i-1)+uj(3,i+1))
       duf(2,i) = dgf(2,i)+dS
       duf(6,i) = dgf(6,i)+dS
    enddo

    call boundary__field(uf,nx,bc)
    call boundary__field(duf,nx,bc)


  end subroutine field__cip


  subroutine ele_cur_cip(uj,up,np,nx,nsp,np2,bc,q,c)

    integer, intent(in)  :: np, nx, nsp, bc
    integer, intent(in)  :: np2(1:nx+bc,nsp)
    real(8), intent(in)  :: q(nsp), c
    real(8), intent(in)  :: up(4,np,1:nx+bc,nsp)
    real(8), intent(out) :: uj(3,-1:nx+2)
    integer :: ii, i, isp
    real(8) :: dx, dxm, gam

    !memory clear
    uj(1:3,0:nx+1) = 0.0D0

    !caluculate erectric current density
    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             gam = 1./dsqrt(1.0+(+up(2,ii,i,isp)*up(2,ii,i,isp) &
                                 +up(3,ii,i,isp)*up(3,ii,i,isp) &
                                 +up(4,ii,i,isp)*up(4,ii,i,isp) &
                                )/(c*c))

             dx = up(1,ii,i,isp)-i
             dxm = 1.-dx

             uj(1  ,i) = uj(1  ,i)+q(isp)*up(2,ii,i,isp)*gam*dxm
             uj(2  ,i) = uj(2  ,i)+q(isp)*up(3,ii,i,isp)*gam*dxm
             uj(3  ,i) = uj(3  ,i)+q(isp)*up(4,ii,i,isp)*gam*dxm

             uj(1,i+1) = uj(1,i+1)+q(isp)*up(2,ii,i,isp)*gam*dx 
             uj(2,i+1) = uj(2,i+1)+q(isp)*up(3,ii,i,isp)*gam*dx 
             uj(3,i+1) = uj(3,i+1)+q(isp)*up(4,ii,i,isp)*gam*dx 
          enddo
       enddo
    enddo

  end subroutine ele_cur_cip


  subroutine particle__solv_cip(gp,up,uf,np,nx,nsp,np2,bc,delt,c,q,r)

    integer, intent(in)  :: np, nx, nsp, bc
    integer, intent(in)  :: np2(1:nx+bc,nsp)
    real(8), intent(in)  :: up(4,np,1:nx+bc,nsp)
    real(8), intent(in)  :: uf(6,0:nx+1)
    real(8), intent(in)  :: delt, c, q(nsp), r(nsp)
    real(8), intent(out) :: gp(4,np,1:nx+bc,nsp)
    integer :: i, ii, isp
    real(8) :: pf(6,np,1:nx+bc)
    real(8) :: uvm(6,np,1:nx+bc)
    real(8) :: dx, dxm
    real(8) :: fac1, fac2, fac2r, fac3, fac3r, gam, txxx, bt2

    do isp=1,nsp

       !interpolate fields to particles
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             dx = up(1,ii,i,isp)-i
             dxm = 1.-dx

             pf(1,ii,i) = +dxm*uf(1,i)+dx*uf(1,i+1)
             pf(2,ii,i) = +dxm*uf(2,i)+dx*uf(2,i+1)
             pf(3,ii,i) = +dxm*uf(3,i)+dx*uf(3,i+1)
             pf(4,ii,i) = +dxm*uf(4,i)+dx*uf(4,i+1)
             pf(5,ii,i) = +dxm*uf(5,i)+dx*uf(5,i+1)
             pf(6,ii,i) = +dxm*uf(6,i)+dx*uf(6,i+1)
          enddo
       enddo

       fac1 = q(isp)/r(isp)*0.5*delt
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             uvm(1,ii,i) = up(2,ii,i,isp)+fac1*pf(4,ii,i)
             uvm(2,ii,i) = up(3,ii,i,isp)+fac1*pf(5,ii,i)
             uvm(3,ii,i) = up(4,ii,i,isp)+fac1*pf(6,ii,i)
          enddo
       enddo

       fac2 = fac1/c
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             gam = dsqrt(1.0+(+uvm(1,ii,i)*uvm(1,ii,i) &
                              +uvm(2,ii,i)*uvm(2,ii,i) &
                              +uvm(3,ii,i)*uvm(3,ii,i))/(c*c))
             fac2r = fac2/gam
             uvm(4,ii,i) = uvm(1,ii,i) &
                  +fac2r*(+uvm(2,ii,i)*pf(3,ii,i)-uvm(3,ii,i)*pf(2,ii,i))
             uvm(5,ii,i) = uvm(2,ii,i) &
                  +fac2r*(+uvm(3,ii,i)*pf(1,ii,i)-uvm(1,ii,i)*pf(3,ii,i))
             uvm(6,ii,i) = uvm(3,ii,i) &
                  +fac2r*(+uvm(1,ii,i)*pf(2,ii,i)-uvm(2,ii,i)*pf(1,ii,i))
          enddo
       enddo
    
       txxx = fac1*fac1
       fac3 = q(isp)*delt/r(isp)
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             bt2 = pf(1,ii,i)*pf(1,ii,i)+pf(2,ii,i)*pf(2,ii,i)+pf(3,ii,i)*pf(3,ii,i)
             gam = dsqrt(1.0+(+uvm(1,ii,i)*uvm(1,ii,i) &
                              +uvm(2,ii,i)*uvm(2,ii,i) &
                              +uvm(3,ii,i)*uvm(3,ii,i))/(c*c))
             fac3r = fac3/(gam+txxx*bt2/gam)
             uvm(1,ii,i) = uvm(1,ii,i) &
                  +fac3r*(+uvm(5,ii,i)*pf(3,ii,i)-uvm(6,ii,i)*pf(2,ii,i))
             uvm(2,ii,i) = uvm(2,ii,i) &
                  +fac3r*(+uvm(6,ii,i)*pf(1,ii,i)-uvm(4,ii,i)*pf(3,ii,i))
             uvm(3,ii,i) = uvm(3,ii,i) &
                  +fac3r*(+uvm(4,ii,i)*pf(2,ii,i)-uvm(5,ii,i)*pf(1,ii,i))
          enddo
       enddo

       do i=1,nx+bc
          do ii=1,np2(i,isp)
             gp(2,ii,i,isp) = uvm(1,ii,i)+fac1*pf(4,ii,i)
             gp(3,ii,i,isp) = uvm(2,ii,i)+fac1*pf(5,ii,i)
             gp(4,ii,i,isp) = uvm(3,ii,i)+fac1*pf(6,ii,i)
          enddo
       enddo
       
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             gam = dsqrt(1.0+(+gp(2,ii,i,isp)*gp(2,ii,i,isp) &
                              +gp(3,ii,i,isp)*gp(3,ii,i,isp) &
                              +gp(4,ii,i,isp)*gp(4,ii,i,isp))/(c*c))
             gp(1,ii,i,isp) = up(1,ii,i,isp)+gp(2,ii,i,isp)*delt/gam
          enddo
       enddo

    enddo

  end subroutine particle__solv_cip
