module particle


  implicit none

  private

  public :: particle__solv


contains


  subroutine particle__solv(up,uf,c,q,r,delt,np,nxgs,nxge,nygs,nyge,nys,nye,nsp,np2)

    integer, intent(in)    :: np, nxgs, nxge, nygs, nyge, nys, nye, nsp
    integer, intent(in)    :: np2(nys:nye,nsp)
    real(8), intent(in)    :: c, q(nsp), r(nsp), delt
    real(8), intent(in)    :: uf(6,nxgs-1:nxge+1,nygs-1:nyge+1)
    real(8), intent(inout) :: up(5,np,nys:nye,nsp)
    integer :: i, j, ii, isp, ih, jh
    real(8) :: dx, dxm, dy, dym
    real(8) :: fac1, fac1r, fac2, fac2r, gam, txxx, bt2
    real(8) :: pf(6)
    real(8) :: uvm(6)

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)

       !interpolate fields to particles
       do j=nys,nye
          do ii=1,np2(j,isp)
             i  = floor(up(1,ii,j,isp))
             ih = floor(up(1,ii,j,isp)-0.5)
             jh = floor(up(2,ii,j,isp)-0.5)

             !Bx at (i+1/2, j)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             pf(1) = +(dxm*uf(1,ih,j  )+dx*uf(1,ih+1,j  ))*dym &
                     +(dxm*uf(1,ih,j+1)+dx*uf(1,ih+1,j+1))*dy

             !By at (i, j+1/2)
             dx = up(1,ii,j,isp)-i
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             pf(2) = +(dxm*uf(2,i,jh  )+dx*uf(2,i+1,jh  ))*dym &
                     +(dxm*uf(2,i,jh+1)+dx*uf(2,i+1,jh+1))*dy

             !Bz at (i, j)
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             pf(3) = +(dxm*uf(3,i,j  )+dx*uf(3,i+1,j  ))*dym &
                     +(dxm*uf(3,i,j+1)+dx*uf(3,i+1,j+1))*dy

             !Ex at (i, j+1/2)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             pf(4) = +(dxm*uf(4,i,jh  )+dx*uf(4,i+1,jh  ))*dym &
                     +(dxm*uf(4,i,jh+1)+dx*uf(4,i+1,jh+1))*dy

             !Ey at (i+1/2, j)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-j
             dym = 1.-dy
             pf(5) = +(dxm*uf(5,ih,j  )+dx*uf(5,ih+1,j  ))*dym &
                     +(dxm*uf(5,ih,j+1)+dx*uf(5,ih+1,j+1))*dy

             !Ez at (i+1/2, j+1/2)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy
             pf(6) = +(dxm*uf(6,ih,jh  )+dx*uf(6,ih+1,jh  ))*dym &
                     +(dxm*uf(6,ih,jh+1)+dx*uf(6,ih+1,jh+1))*dy

             bt2 = pf(1)*pf(1)+pf(2)*pf(2)+pf(3)*pf(3)

             uvm(1) = up(3,ii,j,isp)+fac1*pf(4)
             uvm(2) = up(4,ii,j,isp)+fac1*pf(5)
             uvm(3) = up(5,ii,j,isp)+fac1*pf(6)

             gam = dsqrt(c*c+uvm(1)*uvm(1)+uvm(2)*uvm(2)+uvm(3)*uvm(3))
             fac1r = fac1/gam
             fac2r = fac2/(gam+txxx*bt2/gam)

             uvm(4) = uvm(1)+fac1r*(+uvm(2)*pf(3)-uvm(3)*pf(2))
             uvm(5) = uvm(2)+fac1r*(+uvm(3)*pf(1)-uvm(1)*pf(3))
             uvm(6) = uvm(3)+fac1r*(+uvm(1)*pf(2)-uvm(2)*pf(1))

             uvm(1) = uvm(1)+fac2r*(+uvm(5)*pf(3)-uvm(6)*pf(2))
             uvm(2) = uvm(2)+fac2r*(+uvm(6)*pf(1)-uvm(4)*pf(3))
             uvm(3) = uvm(3)+fac2r*(+uvm(4)*pf(2)-uvm(5)*pf(1))

             up(3,ii,j,isp) = uvm(1)+fac1*pf(4)
             up(4,ii,j,isp) = uvm(2)+fac1*pf(5)
             up(5,ii,j,isp) = uvm(3)+fac1*pf(6)
          enddo
       enddo

    enddo

  end subroutine particle__solv


end module particle
