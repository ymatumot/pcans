module wk
  
  implicit none

  private 

  public :: wk_f

contains
  
  subroutine wk_f(uf,nx,dir,file13,file14)

    integer, intent(in)          :: nx
    real(8), intent(in)          :: uf(6,0:nx+1)
    character(len=*), intent(in) :: dir, file13, file14
    integer :: i

    open(13,file=trim(dir)//trim(file13),status='unknown')
    open(14,file=trim(dir)//trim(file14),status='unknown')

    !save data for w-k diagram
    write(13,'(100000e13.4)')(uf(2,i),i=1,nx)
    write(14,'(100000e13.4)')(uf(3,i),i=1,nx)

  end subroutine wk_f


end module wk
