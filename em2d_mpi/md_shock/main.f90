program main

  use const
  use mpi_set
  use init
  use boundary, only : boundary__particle
  use fio, only : fio__output
  use particle
  use field

  implicit none

  integer :: it

!**********************************************************************c
!
!    two-dimensional electromagnetic plasma simulation code
!
!    written by M Hoshino,  ISAS, 1984/09/12
!    revised  1985/03/08  1985/04/05  1997/05/06
!    revised for CANS    (by Y. Matsumoto, STEL)  2004/06/22
!    re-written in F90   (by Y. Matsumoto, STEL)  2008/10/21
!    MPI parallelization (by Y. Matsumoto, STEL)  2009/4/1
!    2-D code            (by Y. Matsumoto, STEL)  2009/6/5
!
!**********************************************************************c

  call init__set_param
  call fio__output(up,uf,np2,it0)

  do it=1,itmax-it0

     call particle__solv(gp,up,uf,np2)
     call field__fdtd_i(uf,up,gp,np2)
     call boundary__particle(up,np2)

     call init__inject

     if(mod(it+it0,intvl1) == 0) call fio__output(up,uf,np2,it+it0)

  enddo

  call MPI_FINALIZE(nerr)


end program main
