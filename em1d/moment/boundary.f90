module boundary

  implicit none

  private

  public :: boundary__den
  public :: boundary__vel


contains


  subroutine boundary__den(den,nx,bc)

    integer, intent(in)    :: nx, bc
    real(8), intent(inout) :: den(0:nx+1)

    if(bc == 0)then
       !periodic condition
       den(1)  = den(1)+den(nx+1)
       den(nx) = den(nx)+den(0)
    else if(bc==-1) then
       !reflective condition
       den(1)     = den(1)    +den(0)
       den(nx+bc) = den(nx+bc)+den(nx+bc+1)
    endif

  end subroutine boundary__den


  subroutine boundary__vel(vel,nx,bc)

    integer, intent(in)    :: nx, bc
    real(8), intent(inout) :: vel(0:nx+1,3)

    if(bc == 0)then
       !periodic condition
       vel(1 ,1:3) = vel(1 ,1:3)+vel(nx+1,1:3)
       vel(nx,1:3) = vel(nx,1:3)+vel(0   ,1:3)
    else if(bc == -1) then
       !reflective condition
       vel(1    ,1:3) = vel(1    ,1:3)+vel(0      ,1:3)
       vel(nx+bc,1:3) = vel(nx+bc,1:3)+vel(nx+bc+1,1:3)
    endif

  end subroutine boundary__vel


end module boundary
