program main

  use boundary
  use fio
  use particle

  implicit none

  logical           :: lflag=.true.
  integer           :: nproc, ndata, idata, irank
  character(len=64) :: dir
  character(len=64) :: ifile
  real(8)           :: x0, y0, dx, dy
  character(len=64) :: xpos, ypos

  ndata = iargc()
  call getarg(1,dir)
  call getarg(2,xpos)
  call getarg(3,ypos)
  read(xpos,*)x0
  read(ypos,*)y0

  write(*,*)'No. of processes?'
  read(*,*)nproc

  dx = 2.0 !sampling area in the x direction
  dy = 2.0 !sampling area in the y direction

  do idata=4,ndata,nproc
     do irank=0,nproc-1

        call getarg(idata+irank,ifile)
        write(*,'(a)')'reading.....  '//trim(dir)//trim(ifile)

        call fio__input(nproc,dir,ifile)

        call particle__solv(up,uf,c,q,r,0.5*delt,np,nxgs,nxge,nygs,nyge,nys,nye,nsp,np2)
        call boundary__particle(up,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2,bc)
        call fio__psd(up,x0,y0,dx*delx,dy*delx,np,nys,nye,nsp,np2,it0,trim(dir)//'../psd/')

        deallocate(np2)
        deallocate(up)
     enddo
  enddo

end program main
