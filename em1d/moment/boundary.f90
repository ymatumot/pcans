module boundary

  implicit none

  private

  public :: boundary__particle
  public :: boundary__den
  public :: boundary__vel


contains


  subroutine boundary__particle(up,np,nx,nsp,np2,bc)

    integer, intent(in)    :: np, nx, nsp, bc
    integer, intent(inout) :: np2(1:nx+bc,nsp)
    real(8), intent(inout) :: up(4,np,1:nx+bc,nsp)
    integer :: i, ii, iii, isp, ipos
    integer :: cnt(nx+bc,nsp), flag(np,1:nx+bc,nsp), cnt_tmp

    cnt(1:nx+bc,1:nsp) = 0
    flag(1:np,1:nx+bc,1:nsp) = 0

    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)

             ipos = floor(up(1,ii,i,isp))

             if(ipos /= i)then
                if(bc == 0)then
                   if(ipos <= 0)then
                      ipos = nx+ipos
                      up(1,ii,i,isp) = up(1,ii,i,isp)+nx
                   endif
                   if(ipos >= nx+1)then
                      ipos = ipos-nx 
                      up(1,ii,i,isp) = up(1,ii,i,isp)-nx
                   endif
                else if(bc == -1)then
                   if(ipos <= 0)then
                      ipos = 1-ipos
                      up(1,ii,i,isp) = 2.0-up(1,ii,i,isp)
                      up(2,ii,i,isp) = -up(2,ii,i,isp)
                   endif
                   if(ipos >= nx)then
                      ipos = 2*nx-ipos-1 
                      up(1,ii,i,isp) = 2.0*nx-up(1,ii,i,isp)
                      up(2,ii,i,isp) = -up(2,ii,i,isp)
                   endif
                else
                   write(*,*)'choose bc=0 (periodic) or bc=-1 (bounded)'
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
