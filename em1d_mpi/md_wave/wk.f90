module wk
  
  implicit none

  private 

  public :: wk__fio


contains

  
  subroutine wk__fio(uf,nxgs,nxge,nxs,nxe,nxs1,nxe1,dir,file13,file14, &
                     nroot,nrank,mnpr,nrecvc2,ndisp2,ncomw,nerr)

    integer, intent(in)          :: nxgs, nxge, nxs, nxe, nxs1, nxe1
    integer, intent(in)          :: nroot, nrank, mnpr, nrecvc2(:), ndisp2(:), ncomw, nerr
    real(8), intent(in)          :: uf(6,nxs1:nxe1)
    character(len=*), intent(in) :: dir, file13, file14
    integer       :: i
    integer, save :: iflag
    real(8)       :: ufg(6,nxgs:nxge)

    if(iflag /= 1 .and. nrank == nroot)then
       open(13,file=trim(dir)//trim(file13),status='unknown')
       open(14,file=trim(dir)//trim(file14),status='unknown')
       iflag = 1
    endif

    call MPI_GATHERV(uf(1,nxs),6*(nxe-nxs+1),mnpr, &  ! send
                     ufg(1,nxgs),nrecvc2,ndisp2,mnpr,  &  ! receive
                     nroot,ncomw,nerr)

    !save data for w-k diagram
    if(nrank == nroot)then
       write(13,'(10000e13.4)')(ufg(4,i),i=nxgs,nxge)
       write(14,'(10000e13.4)')(ufg(5,i),i=nxgs,nxge)
    endif

  end subroutine wk__fio


end module wk
