module wk
  
  implicit none

  private 

  public :: wk_f

contains
  
  subroutine wk_f(uf,nx,it,dir,file13,file14)

    integer, intent(in)          :: nx,it
    real(8), intent(in)          :: uf(6,0:nx+1)
    character(len=*), intent(in) :: dir, file13, file14
    integer :: i

    if(it == 0)then
       open(13,file=trim(dir)//trim(file13),status='replace')
       open(14,file=trim(dir)//trim(file14),status='replace')
    else
       open(13,file=trim(dir)//trim(file13),position='append')
       open(14,file=trim(dir)//trim(file14),position='append')
    endif

    !save data for w-k diagram
    write(13,'(100000e13.4)')(uf(2,i),i=1,nx)
    write(14,'(100000e13.4)')(uf(3,i),i=1,nx)

  end subroutine wk_f


end module wk
