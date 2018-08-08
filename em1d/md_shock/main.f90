program main

  use const
  use init
  use boundary
  use fio
  use particle
  use field
  use mom_calc

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

  !INITIALIZATIONS
  call init__set_param

  !RECORDING INITIAL DATA
  call fio__output(up,uf,np2,0,file10)
  call fio__energy(up,uf,np2,0,file12)

  do it=1,itmax-it0

     if(mod(it+it0,10) == 0) call fio__progress_bar(it+it0,itmax)

     call particle__solv(gp,up,uf,np2)
     call field__fdtd_i(uf,up,gp,np2)
     call boundary__particle(up,np2)

     call init__inject

     if(mod(it+it0,intvl1) == 0) &
          call fio__output(up,uf,np2,it+it0,file10)
     if(mod(it+it0,intvl2) == 0) &
          call fio__energy(up,uf,np2,it+it0,file12)
     if(mod(it+it0,intvl3) == 0)then
        call mom_calc__accel(gp,up,uf,np2)
        call mom_calc__nvt(den,vel,temp,gp,np2)
        call boundary__mom(den,vel,temp)
        call fio__mom(den,vel,temp,uf,it+it0)
     endif
     if(mod(it+it0,intvl4) == 0) &
          call fio__psd(up,np2,it+it0)

  enddo

end program main
