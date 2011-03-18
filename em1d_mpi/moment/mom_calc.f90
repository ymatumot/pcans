module mom_calc

  implicit none

  private

  public :: mom_calc__den, mom_calc__vel, mom_calc__temp

contains
  
  subroutine mom_calc__den(den,ux,np,nxgs,nxge,nsp,np2,bc)

    integer, intent(in)  :: np, nxgs, nxge, nsp, bc
    integer, intent(in)  :: np2(nxgs:nxge+bc,nsp)
    real(8), intent(in)  :: ux(1,np,nxgs:nxge+bc,nsp)
    real(8), intent(out) :: den(nxgs-1:nxge+1,nsp)
    integer :: ii, i, ih, isp
    real(8) :: dx, dxm, an

    !memory clear
    den(nxgs-1:nxge+1,1:nsp)=0.0D0

    !caluculate number density at x=i+1/2
    do isp=1,nsp
       do i=nxgs,nxge+bc
          do ii=1,np2(i,isp)
             ih = floor(ux(1,ii,i,isp)+0.5)
             dx = ux(1,ii,i,isp)+0.5-ih
             dxm = 1.-dx

             den(ih-1,isp) = den(ih-1,isp)+dxm
             den(ih  ,isp) = den(ih  ,isp)+dx 
          enddo
       enddo
    enddo

  end subroutine mom_calc__den


  subroutine mom_calc__vel(vel,up,c,np,nxgs,nxge,nsp,np2,bc)

    integer, intent(in)  :: np, nxgs, nxge, nsp, bc
    integer, intent(in)  :: np2(nxgs:nxge+bc,nsp)
    real(8), intent(in)  :: c
    real(8), intent(in)  :: up(4,np,nxgs:nxge+bc,nsp)
    real(8), intent(out) :: vel(nxgs-1:nxge+1,3,nsp)
    integer :: ii, i, ih, isp
    real(8) :: dx, dxm, gam

    !memory clear
    vel(nxgs-1:nxge+1,1:3,1:nsp) = 0.0D0

    !caluculate velocity at x=i+1/2
    do isp=1,nsp
       do i=nxgs,nxge+bc
          do ii=1,np2(i,isp)
             gam = 1./dsqrt(1.0+(+up(2,ii,i,isp)*up(2,ii,i,isp) &
                                 +up(3,ii,i,isp)*up(3,ii,i,isp) &
                                 +up(4,ii,i,isp)*up(4,ii,i,isp) &
                                )/(c*c))

             ih = floor(up(1,ii,i,isp)+0.5)
             dx = up(1,ii,i,isp)+0.5-ih
             dxm = 1.-dx

             vel(ih-1,1,isp) = vel(ih-1,1,isp)+up(2,ii,i,isp)*gam*dxm
             vel(ih-1,2,isp) = vel(ih-1,2,isp)+up(3,ii,i,isp)*gam*dxm
             vel(ih-1,3,isp) = vel(ih-1,3,isp)+up(4,ii,i,isp)*gam*dxm
             vel(ih  ,1,isp) = vel(ih  ,1,isp)+up(2,ii,i,isp)*gam*dx 
             vel(ih  ,2,isp) = vel(ih  ,2,isp)+up(3,ii,i,isp)*gam*dx 
             vel(ih  ,3,isp) = vel(ih  ,3,isp)+up(4,ii,i,isp)*gam*dx 
          enddo
       enddo
    enddo

  end subroutine mom_calc__vel


  subroutine mom_calc__temp(temp,up,c,np,nxgs,nxge,nsp,np2,bc)

    integer, intent(in)  :: np, nxgs, nxge, nsp, bc
    integer, intent(in)  :: np2(nxgs:nxge+bc,nsp)
    real(8), intent(in)  :: c
    real(8), intent(in)  :: up(4,np,nxgs:nxge+bc,nsp)
    real(8), intent(out) :: temp(nxgs-1:nxge+1,3,nsp)
    integer :: ii, i, ih, isp
    real(8) :: dx, dxm, gam

    !memory clear
    temp(nxgs-1:nxge+1,1:3,1:nsp) = 0.0D0

    !caluculate velocity at i
    do isp=1,nsp
       do i=nxgs,nxge+bc
          do ii=1,np2(i,isp)
             gam = 1./dsqrt(1.0+(+up(2,ii,i,isp)*up(2,ii,i,isp) &
                                 +up(3,ii,i,isp)*up(3,ii,i,isp) &
                                 +up(4,ii,i,isp)*up(4,ii,i,isp) &
                                )/(c*c))

             ih = floor(up(1,ii,i,isp)+0.5)
             dx = up(1,ii,i,isp)+0.5-ih
             dxm = 1.-dx

             temp(ih-1,1,isp) = temp(ih-1,1,isp)+up(2,ii,i,isp)*up(2,ii,i,isp)*gam*dxm
             temp(ih-1,2,isp) = temp(ih-1,2,isp)+up(3,ii,i,isp)*up(3,ii,i,isp)*gam*dxm
             temp(ih-1,3,isp) = temp(ih-1,3,isp)+up(4,ii,i,isp)*up(4,ii,i,isp)*gam*dxm
             temp(ih  ,1,isp) = temp(ih  ,1,isp)+up(2,ii,i,isp)*up(2,ii,i,isp)*gam*dx 
             temp(ih  ,2,isp) = temp(ih  ,2,isp)+up(3,ii,i,isp)*up(3,ii,i,isp)*gam*dx 
             temp(ih  ,3,isp) = temp(ih  ,3,isp)+up(4,ii,i,isp)*up(4,ii,i,isp)*gam*dx 
          enddo
       enddo
    enddo

  end subroutine mom_calc__temp


end module mom_calc
