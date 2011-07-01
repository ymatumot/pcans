module boundary

  implicit none

  private

  public :: boundary__den
  public :: boundary__vel


contains


  subroutine boundary__den(den,nxgs,nxge,bc)

    integer, intent(in)    :: nxgs, nxge, bc
    real(8), intent(inout) :: den(nxgs-1:nxge+1)

    if(bc == 0)then
       !periodic condition
       den(nxgs) = den(nxgs)+den(nxge+1)
       den(nxge) = den(nxge)+den(nxgs-1)
    else if(bc==-1)then
       !reflective condition
       den(nxgs)    = den(nxgs)   +den(nxgs-1)
       den(nxge+bc) = den(nxge+bc)+den(nxge+bc+1)
    endif

  end subroutine boundary__den


  subroutine boundary__vel(vel,nxgs,nxge,bc)

    integer, intent(in)    :: nxgs, nxge, bc
    real(8), intent(inout) :: vel(nxgs-1:nxge+1,3)

    if(bc == 0)then
       !periodic condition
       vel(nxgs,1:3) = vel(nxgs,1:3)+vel(nxge+1,1:3)
       vel(nxge,1:3) = vel(nxge,1:3)+vel(nxgs-1,1:3)
    else if(bc == -1)then
       !reflective condition
       vel(nxgs   ,1:3) = vel(nxgs   ,1:3)+vel(nxgs-1   ,1:3)
       vel(nxge+bc,1:3) = vel(nxge+bc,1:3)+vel(nxge+bc+1,1:3)
    endif

  end subroutine boundary__vel


end module boundary
