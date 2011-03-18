program main

  use const
  use init
  use boundary
  use fio
  use particle
  use field
  
  implicit none

  integer :: it=0

!**********************************************************************c
!
!    one-dimensional electromagnetic plasma simulation code
!
!    written by M Hoshino,  ISAS, 1984/09/12
!    revised  1985/03/08  1985/04/05  1997/05/06
!    revised for CANS (by Y. Matsumoto, STEL)  2004/06/22
!    re-written in F90 (by Y. Matsumoto, STEL)  2008/10/21
!
!**********************************************************************c

  call init__set_param
  call fio__energy(up,uf,np,nx,nsp,np2,c,r,delt,bc,it,it0,dir,file12)
  call fio__output(up,uf,np,nx,nsp,np2,c,q,r,delt,delx,bc,0,it0,dir,file10)

  do it=1,itmax-it0

     call particle__solv(gp,up,uf,c,q,r,delt,np,nx,nsp,np2,bc)
     call field__fdtd_i(uf,up,gp,np,nx,nsp,np2,bc,q,c,delx,delt,gfac)
     call boundary__particle(up,np,nx,nsp,np2,bc)

     if(mod(it+it0,intvl1) == 0) &
          call fio__output(up,uf,np,nx,nsp,np2,c,q,r,delt,delx,bc,it,it0,dir,file10)
     if(mod(it+it0,intvl2) == 0) &
          call fio__energy(up,uf,np,nx,nsp,np2,c,r,delt,bc,it,it0,dir,file12)

  enddo

end program main
