module boundary

  implicit none

  private

  public :: boundary__particle
  public :: boundary__den
  public :: boundary__vel

contains

  subroutine boundary__particle(up,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2,bc)

    integer, intent(in)    :: np, nys, nye, nxgs, nxge, nygs, nyge, nsp, bc
    integer, intent(inout) :: np2(nys:nye,nsp)
    real(8), intent(inout) :: up(5,np,nys:nye,nsp)
    integer :: j, ii, isp, ipos, jpos

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)

             ipos = floor(up(1,ii,j,isp))
             jpos = floor(up(2,ii,j,isp))

             if(bc==0)then
                if(ipos <= nxgs-1)then
                   up(1,ii,j,isp) = up(1,ii,j,isp)+(nxge-nxgs+1)
                endif
                if(ipos >= nxge+1)then
                   up(1,ii,j,isp) = up(1,ii,j,isp)-(nxge-nxgs+1)
                endif
             else if(bc==-1)then
                if(ipos <= nxgs-1)then
                   up(1,ii,j,isp) = 2.0*nxgs-up(1,ii,j,isp)
                   up(3,ii,j,isp) = -up(3,ii,j,isp)
                endif
                if(ipos >= nxge)then
                   up(1,ii,j,isp) = 2.0*nxge-up(1,ii,j,isp)
                   up(3,ii,j,isp) = -up(3,ii,j,isp)
                endif
             else
                write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                stop
             endif

             if(jpos <= nygs-1)then
                jpos = jpos+(nyge-nygs+1)
                up(2,ii,j,isp) = up(2,ii,j,isp)+(nyge-nygs+1)
             endif
             if(jpos >= nyge+1)then
                jpos = jpos-(nyge-nygs+1)
                up(2,ii,j,isp) = up(2,ii,j,isp)-(nyge-nygs+1)
             endif
          enddo
       enddo
    enddo

  end subroutine boundary__particle


  subroutine boundary__den(den,nxgs,nxge,nygs,nyge,bc)

    integer, intent(in)    :: nxgs, nxge, nygs, nyge, bc
    real(8), intent(inout) :: den(nxgs-1:nxge+1, nygs-1:nyge+1)

    if(bc == 0)then
       !periodic condition
       den(nxgs,nygs-1:nyge+1) = den(nxgs,nygs-1:nyge+1)+den(nxge+1,nygs-1:nyge+1)
       den(nxge,nygs-1:nyge+1) = den(nxge,nygs-1:nyge+1)+den(nxgs-1,nygs-1:nyge+1)
    else if(bc==-1)then
       !reflective condition
       den(nxgs   ,nygs-1:nyge+1) = den(nxgs   ,nygs-1:nyge+1)+den(nxgs-1   ,nygs-1:nyge+1)
       den(nxge+bc,nygs-1:nyge+1) = den(nxge+bc,nygs-1:nyge+1)+den(nxge+bc+1,nygs-1:nyge+1)
    endif

    den(nxgs:nxge+bc,nygs) = den(nxgs:nxge+bc,nygs)+den(nxgs:nxge+bc,nyge+1)
    den(nxgs:nxge+bc,nyge) = den(nxgs:nxge+bc,nyge)+den(nxgs:nxge+bc,nygs-1)

  end subroutine boundary__den


  subroutine boundary__vel(vel,nxgs,nxge,nygs,nyge,bc)

    integer, intent(in)    :: nxgs, nxge, nygs, nyge, bc
    real(8), intent(inout) :: vel(nxgs-1:nxge+1,nygs-1:nyge+1,3)

    if(bc == 0)then
       !periodic condition
       vel(nxgs  ,nygs-1:nyge+1,1:3) = vel(nxgs  ,nygs-1:nyge+1,1:3)+vel(nxge+1,nygs-1:nyge+1,1:3)
       vel(nxge  ,nygs-1:nyge+1,1:3) = vel(nxge  ,nygs-1:nyge+1,1:3)+vel(nxgs-1,nygs-1:nyge+1,1:3)
    else if(bc==-1)then
       !reflective condition
       vel(nxgs   ,nygs-1:nyge+1,1:3) = vel(nxgs   ,nygs-1:nyge+1,1:3)+vel(nxgs-1   ,nygs-1:nyge+1,1:3)
       vel(nxge+bc,nygs-1:nyge+1,1:3) = vel(nxge+bc,nygs-1:nyge+1,1:3)+vel(nxge+bc+1,nygs-1:nyge+1,1:3)
    endif
    vel(nxgs:nxge+bc,nygs,1:3) = vel(nxgs:nxge+bc,nygs,1:3)+vel(nxgs:nxge+bc,nyge+1,1:3)
    vel(nxgs:nxge+bc,nyge,1:3) = vel(nxgs:nxge+bc,nyge,1:3)+vel(nxgs:nxge+bc,nygs-1,1:3)
    
  end subroutine boundary__vel


end module boundary
