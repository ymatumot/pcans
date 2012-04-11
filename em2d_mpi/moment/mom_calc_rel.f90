module mom_calc_rel

  implicit none

  private

  public :: mom_calc__N, mom_calc__T

contains
  
  subroutine mom_calc__N(uN,up,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)

    integer, intent(in)    :: np, nys, nye, nxgs, nxge, nygs, nyge, nsp, np2(nys:nye,nsp)
    real(8), intent(in)    :: c
    real(8), intent(in)    :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: uN(nxgs-1:nxge+1,nygs-1:nyge+1,0:3,1:nsp)
    integer       :: ii, ih, j, jh, isp
    real(8) :: dx, dxm, dy, dym, gam

    ! caluculate particle 4-flow at (i+1/2, j+1/2)
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

             ! N0: number density
             uN(ih  ,jh  ,0,isp) = uN(ih  ,jh  ,0,isp)+dxm*dym
             uN(ih+1,jh  ,0,isp) = uN(ih+1,jh  ,0,isp)+dx *dym
             uN(ih  ,jh+1,0,isp) = uN(ih  ,jh+1,0,isp)+dxm*dy
             uN(ih+1,jh+1,0,isp) = uN(ih+1,jh+1,0,isp)+dx *dy

             ! N1: number flux in x
             uN(ih  ,jh  ,1,isp) = uN(ih  ,jh  ,1,isp)+up(3,ii,j,isp)*gam*dxm*dym
             uN(ih+1,jh  ,1,isp) = uN(ih+1,jh  ,1,isp)+up(3,ii,j,isp)*gam*dx *dym
             uN(ih  ,jh+1,1,isp) = uN(ih  ,jh+1,1,isp)+up(3,ii,j,isp)*gam*dxm*dy
             uN(ih+1,jh+1,1,isp) = uN(ih+1,jh+1,1,isp)+up(3,ii,j,isp)*gam*dx *dy

             ! N2: number flux in y
             uN(ih  ,jh  ,2,isp) = uN(ih  ,jh  ,2,isp)+up(4,ii,j,isp)*gam*dxm*dym
             uN(ih+1,jh  ,2,isp) = uN(ih+1,jh  ,2,isp)+up(4,ii,j,isp)*gam*dx *dym
             uN(ih  ,jh+1,2,isp) = uN(ih  ,jh+1,2,isp)+up(4,ii,j,isp)*gam*dxm*dy
             uN(ih+1,jh+1,2,isp) = uN(ih+1,jh+1,2,isp)+up(4,ii,j,isp)*gam*dx *dy

             ! N3: number flux in z
             uN(ih  ,jh  ,3,isp) = uN(ih  ,jh  ,3,isp)+up(5,ii,j,isp)*gam*dxm*dym
             uN(ih+1,jh  ,3,isp) = uN(ih+1,jh  ,3,isp)+up(5,ii,j,isp)*gam*dx *dym
             uN(ih  ,jh+1,3,isp) = uN(ih  ,jh+1,3,isp)+up(5,ii,j,isp)*gam*dxm*dy
             uN(ih+1,jh+1,3,isp) = uN(ih+1,jh+1,3,isp)+up(5,ii,j,isp)*gam*dx *dy
          enddo
       enddo
    enddo

  end subroutine mom_calc__N


  subroutine mom_calc__T(uT,up,c,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)

    integer, intent(in)    :: np, nys, nye, nxgs, nxge, nygs, nyge, nsp, np2(nys:nye,nsp)
    real(8), intent(in)    :: c
    real(8), intent(in)    :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: uT(nxgs-1:nxge+1,nygs-1:nyge+1,1:10,1:nsp)
    integer :: ii, i, j, ih, jh, isp
    real(8) :: dx, dxm, dy, dym, gam, gam2

    !caluculate stress-energy tensor at (i+1/2, j+1/2)
    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             gam = dsqrt(1.0+(+up(3,ii,j,isp)*up(3,ii,j,isp) &
                              +up(4,ii,j,isp)*up(4,ii,j,isp) &
                              +up(5,ii,j,isp)*up(5,ii,j,isp) &
                             )/(c*c))
             gam2 = 1./gam

             ih = floor(up(1,ii,j,isp)-0.5)
             dx = up(1,ii,j,isp)-0.5-ih
             dxm = 1.-dx
             jh = floor(up(2,ii,j,isp)-0.5)
             dy = up(2,ii,j,isp)-0.5-jh
             dym = 1.-dy

             !T00: energy density
             uT(ih  ,jh  ,1,isp) = uT(ih  ,jh  ,1,isp)+gam*dxm*dym
             uT(ih+1,jh  ,1,isp) = uT(ih+1,jh  ,1,isp)+gam*dx *dym
             uT(ih  ,jh+1,1,isp) = uT(ih  ,jh+1,1,isp)+gam*dxm*dy
             uT(ih+1,jh+1,1,isp) = uT(ih+1,jh+1,1,isp)+gam*dx *dy

             !T01: energy flux in x
             uT(ih  ,jh  ,2,isp) = uT(ih  ,jh  ,2,isp)+up(3,ii,j,isp)*dxm*dym
             uT(ih+1,jh  ,2,isp) = uT(ih+1,jh  ,2,isp)+up(3,ii,j,isp)*dx *dym
             uT(ih  ,jh+1,2,isp) = uT(ih  ,jh+1,2,isp)+up(3,ii,j,isp)*dxm*dy
             uT(ih+1,jh+1,2,isp) = uT(ih+1,jh+1,2,isp)+up(3,ii,j,isp)*dx *dy

             !T02: energy flux in y
             uT(ih  ,jh  ,3,isp) = uT(ih  ,jh  ,3,isp)+up(4,ii,j,isp)*dxm*dym
             uT(ih+1,jh  ,3,isp) = uT(ih+1,jh  ,3,isp)+up(4,ii,j,isp)*dx *dym
             uT(ih  ,jh+1,3,isp) = uT(ih  ,jh+1,3,isp)+up(4,ii,j,isp)*dxm*dy
             uT(ih+1,jh+1,3,isp) = uT(ih+1,jh+1,3,isp)+up(4,ii,j,isp)*dx *dy

             !T03: energy flux in z
             uT(ih  ,jh  ,4,isp) = uT(ih  ,jh  ,4,isp)+up(5,ii,j,isp)*dxm*dym
             uT(ih+1,jh  ,4,isp) = uT(ih+1,jh  ,4,isp)+up(5,ii,j,isp)*dx *dym
             uT(ih  ,jh+1,4,isp) = uT(ih  ,jh+1,4,isp)+up(5,ii,j,isp)*dxm*dy
             uT(ih+1,jh+1,4,isp) = uT(ih+1,jh+1,4,isp)+up(5,ii,j,isp)*dx *dy

             !T11: stress-energy tensor xx
             uT(ih  ,jh  ,5,isp) = uT(ih  ,jh  ,5,isp)+up(3,ii,j,isp)**2*gam2*dxm*dym
             uT(ih+1,jh  ,5,isp) = uT(ih+1,jh  ,5,isp)+up(3,ii,j,isp)**2*gam2*dx *dym
             uT(ih  ,jh+1,5,isp) = uT(ih  ,jh+1,5,isp)+up(3,ii,j,isp)**2*gam2*dxm*dy
             uT(ih+1,jh+1,5,isp) = uT(ih+1,jh+1,5,isp)+up(3,ii,j,isp)**2*gam2*dx *dy

             !T12: stress-energy tensor xy
             uT(ih  ,jh  ,6,isp) = uT(ih  ,jh  ,6,isp)+up(3,ii,j,isp)*up(4,ii,j,isp)*gam2*dxm*dym
             uT(ih+1,jh  ,6,isp) = uT(ih+1,jh  ,6,isp)+up(3,ii,j,isp)*up(4,ii,j,isp)*gam2*dx *dym
             uT(ih  ,jh+1,6,isp) = uT(ih  ,jh+1,6,isp)+up(3,ii,j,isp)*up(4,ii,j,isp)*gam2*dxm*dy
             uT(ih+1,jh+1,6,isp) = uT(ih+1,jh+1,6,isp)+up(3,ii,j,isp)*up(4,ii,j,isp)*gam2*dx *dy

             !T13: stress-energy tensor xz
             uT(ih  ,jh  ,7,isp) = uT(ih  ,jh  ,7,isp)+up(3,ii,j,isp)*up(5,ii,j,isp)*gam2*dxm*dym
             uT(ih+1,jh  ,7,isp) = uT(ih+1,jh  ,7,isp)+up(3,ii,j,isp)*up(5,ii,j,isp)*gam2*dx *dym
             uT(ih  ,jh+1,7,isp) = uT(ih  ,jh+1,7,isp)+up(3,ii,j,isp)*up(5,ii,j,isp)*gam2*dxm*dy
             uT(ih+1,jh+1,7,isp) = uT(ih+1,jh+1,7,isp)+up(3,ii,j,isp)*up(5,ii,j,isp)*gam2*dx *dy

             !T22: stress-energy tensor yy
             uT(ih  ,jh  ,8,isp) = uT(ih  ,jh  ,8,isp)+up(4,ii,j,isp)**2*gam2*dxm*dym
             uT(ih+1,jh  ,8,isp) = uT(ih+1,jh  ,8,isp)+up(4,ii,j,isp)**2*gam2*dx *dym
             uT(ih  ,jh+1,8,isp) = uT(ih  ,jh+1,8,isp)+up(4,ii,j,isp)**2*gam2*dxm*dy
             uT(ih+1,jh+1,8,isp) = uT(ih+1,jh+1,8,isp)+up(4,ii,j,isp)**2*gam2*dx *dy

             !T23: stress-energy tensor yz
             uT(ih  ,jh  ,9,isp) = uT(ih  ,jh  ,9,isp)+up(4,ii,j,isp)*up(5,ii,j,isp)*gam2*dxm*dym
             uT(ih+1,jh  ,9,isp) = uT(ih+1,jh  ,9,isp)+up(4,ii,j,isp)*up(5,ii,j,isp)*gam2*dx *dym
             uT(ih  ,jh+1,9,isp) = uT(ih  ,jh+1,9,isp)+up(4,ii,j,isp)*up(5,ii,j,isp)*gam2*dxm*dy
             uT(ih+1,jh+1,9,isp) = uT(ih+1,jh+1,9,isp)+up(4,ii,j,isp)*up(5,ii,j,isp)*gam2*dx *dy

             !T33: stress-energy tensor zz
             uT(ih  ,jh  ,10,isp) = uT(ih  ,jh  ,10,isp)+up(5,ii,j,isp)**2*gam2*dxm*dym
             uT(ih+1,jh  ,10,isp) = uT(ih+1,jh  ,10,isp)+up(5,ii,j,isp)**2*gam2*dx *dym
             uT(ih  ,jh+1,10,isp) = uT(ih  ,jh+1,10,isp)+up(5,ii,j,isp)**2*gam2*dxm*dy
             uT(ih+1,jh+1,10,isp) = uT(ih+1,jh+1,10,isp)+up(5,ii,j,isp)**2*gam2*dx *dy

          enddo
       enddo
    enddo

  end subroutine mom_calc__T


end module mom_calc_rel
