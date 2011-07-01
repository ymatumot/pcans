module boundary

  implicit none

  private

  public :: boundary__den
  public :: boundary__vel


contains


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
