program main


  use const
  use boundary
  use fio
  use particle
  use mom_calc

  implicit none

  integer :: np2(1:nx+bc,nsp), it0
  integer :: ndata, idata, isp
  real(8) :: up(4,np,1:nx+bc,nsp)
  real(8) :: uf(6,0:nx+1)
  real(8) :: den(0:nx+1,nsp), vel(0:nx+1,3,nsp), temp(0:nx+1,3,nsp)
  real(8) :: c, q(nsp), r(nsp), delt, delx
  character(len=64)  :: dir
  character(len=64) :: ifile

  ndata = iargc()
  call getarg(1,dir)
  
  do idata=2,ndata

     call getarg(idata,ifile)
     write(*,'(a)')'reading....'//trim(dir)//trim(ifile)

     call fio__input(up,uf,c,q,r,delt,delx,np2,it0,np,nx,nsp,bc,dir,ifile)
     call particle__solv(up,uf,c,q,r,0.5*delt,np2,np,nx,nsp,bc)
     call boundary__particle(up,np,nx,nsp,np2,bc)

     call mom_calc__den(den,up(1,1:np,1:nx+bc,1:nsp),np,nx,nsp,np2,bc)
     do isp=1,nsp
        call boundary__den(den(0:nx+1,isp),nx,bc)
     enddo

     call mom_calc__vel(vel,up,c,np,nx,nsp,np2,bc)
     do isp=1,nsp
        call boundary__vel(vel(0:nx+1,1:3,isp),nx,bc)
     enddo

     call mom_calc__temp(temp,up,c,np,nx,nsp,np2,bc)
     do isp=1,nsp
        call boundary__den(temp(0:nx+1,1,isp),nx,bc)
        call boundary__den(temp(0:nx+1,2,isp),nx,bc)
        call boundary__den(temp(0:nx+1,3,isp),nx,bc)
     enddo

     call fio__mom(den,vel,temp,uf,nx,nsp,bc,it0,trim(dir)//'../mom/')
     call fio__psd(up,np,nx,nsp,np2,bc,it0,trim(dir)//'../psd/')
  enddo

end program main
