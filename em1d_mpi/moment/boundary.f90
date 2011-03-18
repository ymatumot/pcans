module boundary

  implicit none

  private

  public :: boundary__particle
  public :: boundary__den
  public :: boundary__vel

contains

  subroutine boundary__particle(up,np,nxgs,nxge,nsp,np2,bc)

    integer, intent(in)    :: np, nxgs, nxge, nsp, bc
    integer, intent(inout) :: np2(nxgs:nxge+bc,nsp)
    real(8), intent(inout) :: up(4,np,nxgs:nxge+bc,nsp)
    integer :: i, ii, iii, isp, ipos
    integer :: cnt(nxgs:nxge+bc,nsp), flag(np,nxgs:nxge+bc,nsp), cnt_tmp

    cnt(nxgs:nxge+bc,1:nsp) = 0
    flag(1:np,nxgs:nxge+bc,1:nsp) = 0

    do isp=1,nsp
       do i=nxgs,nxge+bc
          do ii=1,np2(i,isp)

             ipos = floor(up(1,ii,i,isp))

             if(ipos /= i)then

                if(bc == 0)then
                   if(ipos <= nxgs-1)then
                      ipos = ipos+(nxge-nxgs+1)
                      up(1,ii,i,isp) = up(1,ii,i,isp)+(nxge-nxgs+1)
                   endif
                   if(ipos >= nxge+1)then
                      ipos = ipos-(nxge-nxgs+1)
                      up(1,ii,i,isp) = up(1,ii,i,isp)-(nxge-nxgs+1)
                   endif
                else if(bc == -1)then
                   if(ipos <= nxgs-1)then
                      ipos = 2*nxgs-ipos-1
                      up(1,ii,i,isp) = 2.0*nxgs-up(1,ii,i,isp)
                      up(2,ii,i,isp) = -up(2,ii,i,isp)
                   endif
                   if(ipos >= nxge)then
                      ipos = 2*nxge-ipos-1
                      up(1,ii,i,isp) = 2.0*nxge-up(1,ii,i,isp)
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
       do i=nxgs,nxge+bc
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
       do i=nxgs,nxge+bc
          np2(i,isp) = np2(i,isp)+cnt(i,isp)
          if(np2(i,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(i,isp),i,isp
             stop
          endif
       enddo
    enddo


  end subroutine boundary__particle


  subroutine boundary__den(den,nxgs,nxge,bc)

    integer, intent(in)    :: nxgs, nxge, bc
    real(8), intent(inout) :: den(nxgs-1:nxge+1)

    if(bc == 0)then
       !periodic condition
       den(nxgs) = den(nxgs)+den(nxge+1)
       den(nxge) = den(nxge)+den(nxgs-1)
    else if(bc==-1)then
       !reflective condition
       den(nxgs)    = den(nxgs)   +den(nxgs-1)
       den(nxge+bc) = den(nxge+bc)+den(nxge+bc+1)
    endif

  end subroutine boundary__den


  subroutine boundary__vel(vel,nxgs,nxge,bc)

    integer, intent(in)    :: nxgs, nxge, bc
    real(8), intent(inout) :: vel(nxgs-1:nxge+1,3)

    if(bc == 0)then
       !periodic condition
       vel(nxgs,1:3) = vel(nxgs,1:3)+vel(nxge+1,1:3)
       vel(nxge,1:3) = vel(nxge,1:3)+vel(nxgs-1,1:3)
    else if(bc == -1)then
       !reflective condition
       vel(nxgs   ,1:3) = vel(nxgs   ,1:3)+vel(nxgs-1   ,1:3)
       vel(nxge+bc,1:3) = vel(nxge+bc,1:3)+vel(nxge+bc+1,1:3)
    endif

  end subroutine boundary__vel


end module boundary
