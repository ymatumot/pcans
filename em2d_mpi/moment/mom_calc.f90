module mom_calc

  implicit none

  private

  public :: mom_calc__den, mom_calc__vel, mom_calc__temp

contains
  
  subroutine mom_calc__den(den,up,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)

    integer, intent(in)    :: np, nys, nye, nxgs, nxge, nygs, nyge, nsp, np2(nys:nye,nsp)
    real(8), intent(in)    :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nygs-1:nyge+1,1:nsp)
    integer       :: ii, ih, j, jh, isp
    real(8)       :: dx, dxm, dy, dym

    !caluculate number density at (i+1/2, j+1/2)
    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             ih = floor(up(1,ii,j,isp)-0.5)
             jh = floor(up(2,ii,j,isp)-0.5)

             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy

             den(ih  ,jh  ,isp) = den(ih  ,jh  ,isp)+dxm*dym
             den(ih+1,jh  ,isp) = den(ih+1,jh  ,isp)+dx *dym
             den(ih  ,jh+1,isp) = den(ih  ,jh+1,isp)+dxm*dy
             den(ih+1,jh+1,isp) = den(ih+1,jh+1,isp)+dx *dy
          enddo
       enddo
    enddo

  end subroutine mom_calc__den


  subroutine mom_calc__vel(vel,up,c,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)

    integer, intent(in)    :: np, nys, nye, nxgs, nxge, nygs, nyge, nsp, np2(nys:nye,nsp)
    real(8), intent(in)    :: c
    real(8), intent(in)    :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: vel(nxgs-1:nxge+1,nygs-1:nyge+1,1:3,1:nsp)
    integer :: ii, i, j, ih, jh, isp
    real(8) :: dx, dxm, dy, dym, gam

    !caluculate velocity at (i+1/2, j+1/2)
    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             gam = 1./dsqrt(1.0+(+up(3,ii,j,isp)*up(3,ii,j,isp) &
                                 +up(4,ii,j,isp)*up(4,ii,j,isp) &
                                 +up(5,ii,j,isp)*up(5,ii,j,isp) &
                                )/(c*c))

             ih = floor(up(1,ii,j,isp)-0.5)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             jh = floor(up(2,ii,j,isp)-0.5)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy

             !Vx
             vel(ih  ,jh  ,1,isp) = vel(ih  ,jh  ,1,isp)+up(3,ii,j,isp)*gam*dxm*dym
             vel(ih+1,jh  ,1,isp) = vel(ih+1,jh  ,1,isp)+up(3,ii,j,isp)*gam*dx *dym
             vel(ih  ,jh+1,1,isp) = vel(ih  ,jh+1,1,isp)+up(3,ii,j,isp)*gam*dxm*dy
             vel(ih+1,jh+1,1,isp) = vel(ih+1,jh+1,1,isp)+up(3,ii,j,isp)*gam*dx *dy

             !Vy
             vel(ih  ,jh  ,2,isp) = vel(ih  ,jh  ,2,isp)+up(4,ii,j,isp)*gam*dxm*dym
             vel(ih+1,jh  ,2,isp) = vel(ih+1,jh  ,2,isp)+up(4,ii,j,isp)*gam*dx *dym
             vel(ih  ,jh+1,2,isp) = vel(ih  ,jh+1,2,isp)+up(4,ii,j,isp)*gam*dxm*dy
             vel(ih+1,jh+1,2,isp) = vel(ih+1,jh+1,2,isp)+up(4,ii,j,isp)*gam*dx *dy

             !Vz
             vel(ih  ,jh  ,3,isp) = vel(ih  ,jh  ,3,isp)+up(5,ii,j,isp)*gam*dxm*dym
             vel(ih+1,jh  ,3,isp) = vel(ih+1,jh  ,3,isp)+up(5,ii,j,isp)*gam*dx *dym
             vel(ih  ,jh+1,3,isp) = vel(ih  ,jh+1,3,isp)+up(5,ii,j,isp)*gam*dxm*dy
             vel(ih+1,jh+1,3,isp) = vel(ih+1,jh+1,3,isp)+up(5,ii,j,isp)*gam*dx *dy
          enddo
       enddo
    enddo

  end subroutine mom_calc__vel


  subroutine mom_calc__temp(temp,up,c,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)

    integer, intent(in)    :: np, nys, nye, nxgs, nxge, nygs, nyge, nsp, np2(nys:nye,nsp)
    real(8), intent(in)    :: c
    real(8), intent(in)    :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: temp(nxgs-1:nxge+1,nygs-1:nyge+1,1:3,1:nsp)
    integer :: ii, i, j, ih, jh, isp
    real(8) :: dx, dxm, dy, dym, gam

    !caluculate velocity at (i+1/2, j+1/2)
    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             gam = 1./dsqrt(1.0+(+up(3,ii,j,isp)*up(3,ii,j,isp) &
                                 +up(4,ii,j,isp)*up(4,ii,j,isp) &
                                 +up(5,ii,j,isp)*up(5,ii,j,isp) &
                                )/(c*c))

             ih = floor(up(1,ii,j,isp)-0.5)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             jh = floor(up(2,ii,j,isp)-0.5)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy

             !Txx
             temp(ih  ,jh  ,1,isp) = temp(ih  ,jh  ,1,isp)+up(3,ii,j,isp)**2*gam*dxm*dym
             temp(ih+1,jh  ,1,isp) = temp(ih+1,jh  ,1,isp)+up(3,ii,j,isp)**2*gam*dx *dym
             temp(ih  ,jh+1,1,isp) = temp(ih  ,jh+1,1,isp)+up(3,ii,j,isp)**2*gam*dxm*dy
             temp(ih+1,jh+1,1,isp) = temp(ih+1,jh+1,1,isp)+up(3,ii,j,isp)**2*gam*dx *dy

             !Tyy
             temp(ih  ,jh  ,2,isp) = temp(ih  ,jh  ,2,isp)+up(4,ii,j,isp)**2*gam*dxm*dym
             temp(ih+1,jh  ,2,isp) = temp(ih+1,jh  ,2,isp)+up(4,ii,j,isp)**2*gam*dx *dym
             temp(ih  ,jh+1,2,isp) = temp(ih  ,jh+1,2,isp)+up(4,ii,j,isp)**2*gam*dxm*dy
             temp(ih+1,jh+1,2,isp) = temp(ih+1,jh+1,2,isp)+up(4,ii,j,isp)**2*gam*dx *dy

             !Tzz
             temp(ih  ,jh  ,3,isp) = temp(ih  ,jh  ,3,isp)+up(5,ii,j,isp)**2*gam*dxm*dym
             temp(ih+1,jh  ,3,isp) = temp(ih+1,jh  ,3,isp)+up(5,ii,j,isp)**2*gam*dx *dym
             temp(ih  ,jh+1,3,isp) = temp(ih  ,jh+1,3,isp)+up(5,ii,j,isp)**2*gam*dxm*dy
             temp(ih+1,jh+1,3,isp) = temp(ih+1,jh+1,3,isp)+up(5,ii,j,isp)**2*gam*dx *dy
          enddo
       enddo
    enddo

  end subroutine mom_calc__temp


end module mom_calc
