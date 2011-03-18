!---- Charge Conservation Method for Jx, Jy ----
  subroutine ele_cur2(uj,up,gp, &
                      np,nsp,np2,nxs,nxe,nys,nye,bc,q,c,delt)

    integer, intent(in)  :: np, nsp, nxs, nxe, nys, nye, bc
    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: q(nsp), c, delt
    real(8), intent(in)  :: up(5,np,nys:nye,nsp), gp(5,np,nys:nye,nsp)
    real(8), intent(out) :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    
    integer :: ii, j, isp
    integer :: i1 ,i2 ,j1 ,j2, ih, jh
    real(8) :: x1, x2, y1, y2, xr, yr, qvx1, qvx2, qvy1, qvy2, idelt
    real(8) :: dx1, dx2, dy1, dy2, dxm1, dxm2, dym1, dym2, dx, dxm, dy, dym

    idelt = 1.D0/delt

    uj(1:3,nxs-2:nxe+2,nys-2:nye+2) = 0.D0

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)

             x1  = up(1,ii,j,isp)
             x2  = gp(1,ii,j,isp)
             y1  = up(2,ii,j,isp)
             y2  = gp(2,ii,j,isp)
             i1  = floor(x1-0.5)
             j1  = floor(y1-0.5)
             i2  = floor(x2-0.5)
             j2  = floor(y2-0.5)

             if(i1==i2) then
                xr = 0.5*(x1+x2)
             else
                xr = max(i1,i2)*delx
             endif
             if(j1==j2) then
                yr = 0.5*(y1+y2)
             else
                yr = max(j1,j2)*delx
             endif

             qvx1 = q(isp)*(xr-x1)*idelt
             qvx2 = q(isp)*up(3,ii,j,isp)-qvx1
             qvy1 = q(isp)*(yr-y1)*idelt
             qvy2 = q(isp)*up(4,ii,j,isp)-qvy1

             dx1  = 0.5*(x1+xr)-0.5-i1
             dx2  = 0.5*(xr+x2)-0.5-i2
             dy1  = 0.5*(y1+yr)-0.5-j1
             dy2  = 0.5*(yr+y2)-0.5-j2
             dxm1 = dx1-1.
             dxm2 = dx2-1.
             dym1 = dy1-1.
             dym2 = dy2-1.

             uj(1,i1  ,j1) = uj(1,i1  ,j1)+qvx1*dym1
             uj(1,i1,j1+1) = uj(1,i1,j1+1)+qvx1*dy1
             uj(1,i2  ,j2) = uj(1,i2  ,j2)+qvx2*dym2
             uj(1,i2,j2+1) = uj(1,i2,j2+1)+qvx2*dy2
             uj(2,i1  ,j1) = uj(2,i1  ,j1)+qvy1*dxm1
             uj(2,i1+1,j1) = uj(2,i1+1,j1)+qvy1*dx1
             uj(2,i2  ,j2) = uj(2,i2  ,j2)+qvy2*dxm2
             uj(2,i2+1,j2) = uj(2,i2+1,j2)+qvy2*dx2

             dx  = x1-0.5-i1
             dxm = 1.-dx
             dy  = y1-0.5-j1
             dym = 1.-dy
             uj(3,i1    ,j1) = uj(3,i1    ,j1)+q(isp)*up(5,ii,j,isp)*dxm*dym
             uj(3,i1+1  ,j1) = uj(3,i1+1  ,j1)+q(isp)*up(5,ii,j,isp)*dy *dym
             uj(3,i1  ,j1+1) = uj(3,i1  ,j1+1)+q(isp)*up(5,ii,j,isp)*dxm*dy
             uj(3,i1+1,j1+1) = uj(3,i1+1,j1+1)+q(isp)*up(5,ii,j,isp)*dy *dy
          enddo
       enddo
    enddo

  end subroutine ele_cur2
