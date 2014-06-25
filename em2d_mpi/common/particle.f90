module particle

  implicit none

  private

  public :: particle__solv


contains


  subroutine particle__solv(gp,up,uf,                        &
                            np,nsp,np2,nxs,nxe,nys,nye,nsfo, &
                            c,q,r,delt,delx)
                            

    use shape_function, only : sf

    integer, intent(in)  :: np, nxs, nxe, nys, nye, nsp, nsfo
    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp)
    real(8), intent(in)  :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    real(8), intent(in)  :: c, q(nsp), r(nsp), delt, delx
    real(8), intent(out) :: gp(5,np,nys:nye,nsp)
    integer :: i, j, ii, isp, ih, ip, jp
    real(8) :: idelx, fac1, fac1r, fac2, fac2r, gam, txxx, bt2
    real(8) :: sh(-2:2,2)
    real(8) :: pf(6)
    real(8) :: uvm(6)
    real(8) :: tmp(1:6,nxs-1:nxe+1,nys-1:nye+1)

    idelx = 1./delx

    !fields at (i+1/2, j+1/2)
    do j=nys-1,nye+1
    do i=nxs-1,nxe+1
       tmp(1,i,j) = 0.5*(+uf(1,i,j)+uf(1,i,j+1))
       tmp(2,i,j) = 0.5*(+uf(2,i,j)+uf(2,i+1,j))
       tmp(3,i,j) = 0.25*(+uf(3,i,j)  +uf(3,i+1,j) &
                          +uf(3,i,j+1)+uf(3,i+1,j+1))
       tmp(4,i,j) = 0.5*(+uf(4,i,j)+uf(4,i+1,j))
       tmp(5,i,j) = 0.5*(+uf(5,i,j)+uf(5,i,j+1))
       tmp(6,i,j) = uf(6,i,j)
    enddo
    enddo

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)
       do j=nys,nye
          do ii=1,np2(j,isp)

             pf(1:6) = 0.D0

             ih = floor(up(1,ii,j,isp))

             sh(-2:2,1) = sf(ih,up(1,ii,j,isp)*idelx-0.5,nsfo)
             sh(-2:2,2) = sf(j ,up(2,ii,j,isp)*idelx-0.5,nsfo)

             !Bx at (i+1/2, j)
             do jp=-1,1
             do ip=-1,1
                pf(1) = pf(1)+tmp(1,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !By at (i, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(2) = pf(2)+tmp(2,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Bz at (i, j)
             do jp=-1,1
             do ip=-1,1
                pf(3) = pf(3)+tmp(3,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ex at (i, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(4) = pf(4)+tmp(4,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ey at (i+1/2, j)
             do jp=-1,1
             do ip=-1,1
                pf(5) = pf(5)+tmp(5,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ez at (i+1/2, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(6) = pf(6)+tmp(6,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

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

             gp(3,ii,j,isp) = uvm(1)+fac1*pf(4)
             gp(4,ii,j,isp) = uvm(2)+fac1*pf(5)
             gp(5,ii,j,isp) = uvm(3)+fac1*pf(6)

             gam = 1./dsqrt(1.0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                 +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                 +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c))
             gp(1,ii,j,isp) = up(1,ii,j,isp)+gp(3,ii,j,isp)*delt*gam
             gp(2,ii,j,isp) = up(2,ii,j,isp)+gp(4,ii,j,isp)*delt*gam
          enddo
       enddo
    enddo

  end subroutine particle__solv


end module particle
