module boundary

  implicit none

  private

  public :: boundary__init
  public :: boundary__field
  public :: boundary__particle
  public :: boundary__curre
  public :: boundary__phi
  public :: boundary__mom

  logical, save :: is_init = .false.
  integer, save :: np, nx, nsp, bc
  real(8), save :: d_delx, delx


contains


  subroutine boundary__init(npin,nxin,nspin,bcin,delxin)

    integer, intent(in) :: npin, nxin, nspin, bcin
    real(8), intent(in) :: delxin
    
    np  = npin
    nx  = nxin
    nsp = nspin
    bc  = bcin
    delx = delxin
    d_delx = 1./delxin
    is_init = .true.

  end subroutine boundary__init


  subroutine boundary__particle(up,np2)

    integer, intent(inout) :: np2(1:nx+bc,nsp)
    real(8), intent(inout) :: up(4,np,1:nx+bc,nsp)
    logical, save          :: lflag=.true.
    integer :: i, ii, iii, isp, ipos, cnt_tmp
    integer :: cnt(1:nx+bc,nsp), flag(np,1:nx+bc,nsp)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary__init()'
       stop
    endif

    cnt(1:nx+bc,1:nsp) = 0
    flag(1:np,1:nx+bc,1:nsp) = 0
    
    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)

             ipos = floor(up(1,ii,i,isp)*d_delx)

             if(ipos /= i)then
                if(bc == 0)then
                   if(ipos <= 0)then
                      ipos = nx+ipos
                      up(1,ii,i,isp) = up(1,ii,i,isp)+nx*delx
                   endif
                   if(ipos >= nx+1)then
                      ipos = ipos-nx 
                      up(1,ii,i,isp) = up(1,ii,i,isp)-nx*delx
                   endif
                else if(bc == -1)then
                   if(ipos <= 0)then
                      ipos = 1-ipos
                      up(1,ii,i,isp) = 2.*delx-up(1,ii,i,isp)
                      up(2,ii,i,isp) = -up(2,ii,i,isp)
                   endif
                   if(ipos >= nx)then
                      ipos = 2*nx-ipos-1
                      up(1,ii,i,isp) = 2.0*nx*delx-up(1,ii,i,isp)
                      up(2,ii,i,isp) = -up(2,ii,i,isp)
                   endif
                else
                   write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                   stop
                endif
                cnt(ipos,isp) = cnt(ipos,isp)+1
                up(1:4,np2(ipos,isp)+cnt(ipos,isp),ipos,isp) = up(1:4,ii,i,isp)
                flag(ii,i,isp) = 1
             endif

          enddo
       enddo
    enddo

    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             if(flag(ii,i,isp) == 1)then
                cnt_tmp = cnt(i,isp)
                loop1: do iii=np2(i,isp)+cnt_tmp,ii,-1
                   cnt(i,isp) = cnt(i,isp)-1
                   if(flag(iii,i,isp) /= 1)then
                      up(1:4,ii,i,isp) = up(1:4,iii,i,isp)
                      exit loop1
                   endif
                enddo loop1
             endif
          enddo
       enddo
    enddo

    do isp=1,nsp
       do i=1,nx+bc
          np2(i,isp) = np2(i,isp)+cnt(i,isp)
          if(np2(i,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(i,isp),i,isp
             stop
          endif
       enddo
    enddo
          
  end subroutine boundary__particle


  subroutine boundary__field(uf)

    real(8), intent(inout) :: uf(6,0:nx+1)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary__init()'
       stop
    endif

    if(bc == 0)then
       !periodic condition
       uf(1:6,0) = uf(1:6,nx)
       uf(1:6,nx+1) = uf(1:6,1)
    else if(bc == -1)then
       !reflective condition
       uf(1,0)   = +uf(1,1)
       uf(2:3,0) = +uf(2:3,2)
       uf(4,0)   = +uf(4,2)
       uf(5:6,0) = -uf(5:6,1)

       uf(1,nx)     = +uf(1,nx-1)
       uf(2:3,nx+1) = +uf(2:3,nx-1)
       uf(4,nx+1)   = +uf(4,nx-1)
       uf(5:6,nx)   = -uf(5:6,nx-1)
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__field


  subroutine boundary__curre(uj)

    real(8), intent(inout) :: uj(3,-1:nx+2)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary__init()'
       stop
    endif

    !PERIODIC CONDITION
    if(bc == 0)then
       uj(1:3,1)    = uj(1:3,1)+uj(1:3,nx+1)
       uj(1:3,2)    = uj(1:3,2)+uj(1:3,nx+2)
       uj(1:3,nx-1) = uj(1:3,nx-1)+uj(1:3,-1)
       uj(1:3,nx)   = uj(1:3,nx)+uj(1:3,0)
    !REFLECTIVE CONDITION
    else if(bc == -1)then
       uj(1  ,2)  = uj(1  ,2)+uj(1  ,0)
       uj(2:3,1)  = uj(2:3,1)-uj(2:3,0)
       uj(2:3,2)  = uj(2:3,2)-uj(2:3,-1)

       uj(1  ,nx-1) = uj(1  ,nx-1)+uj(1  ,nx+1)
       uj(2:3,nx-1) = uj(2:3,nx-1)-uj(2:3,nx)
       uj(1  ,nx-2) = uj(1  ,nx-2)+uj(1  ,nx+2)
       uj(2:3,nx-2) = uj(2:3,nx-2)-uj(2:3,nx+1)
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__curre


  subroutine boundary__phi(x)

    real(8), intent(inout) :: x(0:nx+1)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary__init()'
       stop
    endif

    !PERIODIC CONDITION
    if(bc == 0)then
       x(0) = x(nx)
       x(nx+1) = x(1)
    !REFLECTIVE CONDITION
    else if(bc == -1)then
       x(0)  = x(2)
       x(nx+1) = x(nx-1)
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__phi


  subroutine boundary__mom(den,vel,temp)

    real(8), intent(inout) :: den(0:nx+1,nsp)
    real(8), intent(inout) :: vel(0:nx+1,3,nsp)
    real(8), intent(inout) :: temp(0:nx+1,3,nsp)
    integer :: isp

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary__init()'
       stop
    endif

    !PERIODIC CONDITION
    if(bc == 0)then
       do isp=1,nsp
          den(1,isp)  = den(1,isp)+den(nx+1,isp)
          den(nx,isp) = den(nx,isp)+den(0,isp)
          vel(1 ,1:3,isp) = vel(1 ,1:3,isp)+vel(nx+1,1:3,isp)
          vel(nx,1:3,isp) = vel(nx,1:3,isp)+vel(0   ,1:3,isp)
          temp(1 ,1:3,isp) = temp(1 ,1:3,isp)+temp(nx+1,1:3,isp)
          temp(nx,1:3,isp) = temp(nx,1:3,isp)+temp(0   ,1:3,isp)
       enddo
    !REFLECTIVE CONDITION
    else if(bc==-1) then
       do isp=1,nsp
          den(1,isp)     = den(1,isp)    +den(0,isp)
          den(nx+bc,isp) = den(nx+bc,isp)+den(nx+bc+1,isp)
          vel(1    ,1:3,isp) = vel(1    ,1:3,isp)+vel(0      ,1:3,isp)
          vel(nx+bc,1:3,isp) = vel(nx+bc,1:3,isp)+vel(nx+bc+1,1:3,isp)
          temp(1    ,1:3,isp) = temp(1    ,1:3,isp)+temp(0      ,1:3,isp)
          temp(nx+bc,1:3,isp) = temp(nx+bc,1:3,isp)+temp(nx+bc+1,1:3,isp)
       enddo
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__mom


end module boundary
