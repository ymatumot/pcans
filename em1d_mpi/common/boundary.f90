module boundary

  implicit none

  private

  public  :: boundary__field
  public  :: boundary__particle
  public  :: boundary__curre
  public  :: boundary__charge


contains


  subroutine boundary__particle(up,                                            &
                                np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bc,bcp, &
                                nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)
    integer, intent(in)    :: np, nsp, nxgs, nxge, nxs, nxe, nxs1, nxe1, bc, bcp
    integer, intent(in)    :: nup, ndown, mnpi, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    integer, intent(inout) :: np2(nxs:nxe+bcp,nsp)
    real(8), intent(inout) :: up(4,np,nxs:nxe+bcp,nsp)
    integer :: i, ii, iii, isp, ipos 
    integer :: cnt(nxs1:nxe1+bcp,nsp), cnt2(nxs:nxe+bcp,nsp), cnt_tmp
    integer :: flag(np,nxs:nxe+bcp,nsp)
    real(8) :: bff_ptcl(np*4,nxs1:nxe1+bcp,nsp)

    cnt(nxs1:nxe1+bcp,1:nsp) = 0
    cnt2(nxs:nxe+bcp,1:nsp) = 0
    flag(1:np,nxs:nxe+bcp,1:nsp) = 0

    do isp=1,nsp
       do i=nxs,nxe+bcp
          do ii=1,np2(i,isp)

             ipos = floor(up(1,ii,i,isp))

             if(ipos /= i)then
                if(bc == 0)then
                   if(ipos <= nxgs-1)then
                      up(1,ii,i,isp) = up(1,ii,i,isp)+(nxge-nxgs+1)
                   endif
                   if(ipos >= nxge+1)then
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
                   write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                   stop
                endif

                bff_ptcl(1+4*cnt(ipos,isp),ipos,isp) = up(1,ii,i,isp)
                bff_ptcl(2+4*cnt(ipos,isp),ipos,isp) = up(2,ii,i,isp)
                bff_ptcl(3+4*cnt(ipos,isp),ipos,isp) = up(3,ii,i,isp)
                bff_ptcl(4+4*cnt(ipos,isp),ipos,isp) = up(4,ii,i,isp)
                cnt(ipos,isp) = cnt(ipos,isp)+1
                cnt2(i,isp) = cnt2(i,isp)+1
                flag(cnt2(i,isp),i,isp) = ii
             endif

          enddo
       enddo
    enddo

    do isp=1,nsp
       !transfer to rank-1
       cnt_tmp = 0
       call MPI_SENDRECV(cnt(nxs1,isp),1,mnpi,ndown,100, &
                         cnt_tmp      ,1,mnpi,nup  ,100, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1               ,nxs1,isp),4*cnt(nxs1,isp),mnpr,ndown,101, &
                         bff_ptcl(4*cnt(nxe,isp)+1,nxe ,isp),4*cnt_tmp      ,mnpr,nup  ,101, &
                         ncomw,nstat,nerr)
       cnt(nxe+bcp,isp) = cnt(nxe+bcp,isp)+cnt_tmp

       !transfer to rank+1
       cnt_tmp = 0
       call MPI_SENDRECV(cnt(nxe1+bcp,isp),1,mnpi,nup  ,200, &
                         cnt_tmp          ,1,mnpi,ndown,200, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1               ,nxe1+bcp,isp),4*cnt(nxe1+bcp,isp),mnpr,nup  ,201, &
                         bff_ptcl(4*cnt(nxs,isp)+1,nxs     ,isp),4*cnt_tmp          ,mnpr,ndown,201, &
                         ncomw,nstat,nerr)
       cnt(nxs,isp) = cnt(nxs,isp)+cnt_tmp
    enddo

    do isp=1,nsp
       do i=nxs,nxe+bcp
          iii=0
          cnt_tmp = cnt2(i,isp)
          loop1 :do ii=1,cnt2(i,isp)
             if(cnt(i,isp) == 0)then
                if(np2(i,isp) < flag(ii,i,isp)) exit loop1
                do while(np2(i,isp) == flag(cnt_tmp,i,isp))
                   np2(i,isp) = np2(i,isp)-1
                   if(np2(i,isp) < flag(ii,i,isp)) exit loop1
                   cnt_tmp = cnt_tmp-1
                enddo
                up(1:4,flag(ii,i,isp),i,isp) = up(1:4,np2(i,isp),i,isp)
                np2(i,isp) = np2(i,isp)-1
             else
                up(1,flag(ii,i,isp),i,isp) = bff_ptcl(1+4*iii,i,isp)
                up(2,flag(ii,i,isp),i,isp) = bff_ptcl(2+4*iii,i,isp)
                up(3,flag(ii,i,isp),i,isp) = bff_ptcl(3+4*iii,i,isp)
                up(4,flag(ii,i,isp),i,isp) = bff_ptcl(4+4*iii,i,isp)
                iii = iii+1
                cnt(i,isp) = cnt(i,isp)-1
             endif
          enddo loop1
          
          if(cnt(i,isp) > 0)then
             do ii=1,cnt(i,isp)
                up(1,np2(i,isp)+ii,i,isp) = bff_ptcl(4*iii+1+4*(ii-1),i,isp)
                up(2,np2(i,isp)+ii,i,isp) = bff_ptcl(4*iii+2+4*(ii-1),i,isp)
                up(3,np2(i,isp)+ii,i,isp) = bff_ptcl(4*iii+3+4*(ii-1),i,isp)
                up(4,np2(i,isp)+ii,i,isp) = bff_ptcl(4*iii+4+4*(ii-1),i,isp)
             enddo
          endif
       enddo
    enddo

    do isp=1,nsp
       do i=nxs,nxe+bcp
          np2(i,isp) = np2(i,isp)+cnt(i,isp)
          if(np2(i,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(i,isp),i,isp
             stop
          endif
       enddo
    enddo

  end subroutine boundary__particle


  subroutine boundary__field(uf,                       &
                             nxs,nxe,nxs1,nxe1,bc,bcp, &
                             nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nxs1, nxe1, bc, bcp
    integer, intent(in)    :: nup, ndown, nroot, nproc, nrank, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uf(6,nxs1:nxe1)
    integer                :: ieq
    real(8)                :: bff_snd(1:6), bff_rcv(1:6)

    bff_snd(1:6) = uf(1:6,nxs)
    call MPI_SENDRECV(bff_snd(1),6,mnpr,ndown,110, &
                      bff_rcv(1),6,mnpr,nup  ,110, &
                      ncomw,nstat,nerr)
    uf(1  ,nxe1+bcp) = bff_rcv(1)
    uf(2:4,nxe1)     = bff_rcv(2:4)
    uf(5:6,nxe1+bcp) = bff_rcv(5:6)

    bff_snd(1)   = uf(1  ,nxe+bcp)
    bff_snd(2:4) = uf(2:4,nxe)
    bff_snd(5:6) = uf(5:6,nxe+bcp)
    call MPI_SENDRECV(bff_snd(1),6,mnpr,nup  ,100, &
                      bff_rcv(1),6,mnpr,ndown,100, &
                      ncomw,nstat,nerr)
    uf(1:6,nxs1) = bff_rcv(1:6)

    if(bc == -1)then
       !reflective condition
       if(nrank == nroot)then
          uf(1,nxs1)   = +uf(1,nxs)
          uf(2:4,nxs1) = +uf(2:4,nxs+1)
          uf(5:6,nxs1) = -uf(5:6,nxs)
       endif
       if(nrank == nproc-1)then
          uf(1  ,nxe1+bcp) = +uf(1  ,nxe1+bcp-1)
          uf(2:4,nxe1)     = +uf(2:4,nxe-1)
          uf(5:6,nxe1+bcp) = -uf(5:6,nxe1+bcp-1)
       endif
    endif

  end subroutine boundary__field


  subroutine boundary__curre(uj,nxs,nxe,nxs1,nxe1,bc, &
                             nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nxs1, nxe1, bc
    integer, intent(in)    :: nup, ndown, nroot, nproc, nrank, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: uj(3,nxs1-1:nxe1+1)
    integer                :: l
    real(8)                :: jtmp, bff_rcv(6), bff_snd(6)

    bff_rcv(1:6) = 0.D0
    bff_snd(1:3) = uj(1:3,nxs1-1)
    bff_snd(4:6) = uj(1:3,nxs1)
    call MPI_SENDRECV(bff_snd(1),6,mnpr,ndown,110, &
                      bff_rcv(1),6,mnpr,nup  ,110, &
                      ncomw,nstat,nerr)
    uj(1:3,nxe-1) = uj(1:3,nxe-1)+bff_rcv(1:3)
    uj(1:3,nxe)   = uj(1:3,nxe  )+bff_rcv(4:6)

    bff_rcv(1:6) = 0.D0
    bff_snd(1:3) = uj(1:3,nxe1)
    bff_snd(4:6) = uj(1:3,nxe1+1)
    call MPI_SENDRECV(bff_snd(1),6,mnpr,nup  ,100, &
                      bff_rcv(1),6,mnpr,ndown,100, &
                      ncomw,nstat,nerr)
    uj(1:3,nxs)   = uj(1:3,nxs)  +bff_rcv(1:3)
    uj(1:3,nxs+1) = uj(1:3,nxs+1)+bff_rcv(4:6)

    if(bc == -1)then
       !reflective condition
       if(nrank == nroot)then
          uj(2:3,nxs)   = uj(2:3,nxs)  -uj(2:3,nxs1)
          uj(1  ,nxs+1) = uj(1  ,nxs+1)+uj(1  ,nxs1)
          uj(2:3,nxs+1) = uj(2:3,nxs+1)-uj(2:3,nxs1-1)
          uj(1  ,nxs+2) = uj(1  ,nxs+2)+uj(1  ,nxs1-1)
       endif
       if(nrank == nproc-1)then
          uj(1  ,nxe-2) = uj(1  ,nxe-2)+uj(1  ,nxe1+1)
          uj(2:3,nxe-2) = uj(2:3,nxe-2)-uj(2:3,nxe1)
          uj(1  ,nxe-1) = uj(1  ,nxe-1)+uj(1  ,nxe1)
          uj(2:3,nxe-1) = uj(2:3,nxe-1)-uj(2:3,nxe)
       endif
    endif

  end subroutine boundary__curre


  subroutine boundary__charge(cden,nxs,nxe,nxs1,nxe1,bc, &
                              nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

    integer, intent(in)    :: nxs, nxe, nxs1, nxe1, bc
    integer, intent(in)    :: nup, ndown, nroot, nproc, nrank, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(inout) :: cden(nxs1-1:nxe1+1)
    real(8)                :: ctmp(2)


    ctmp(1:2) = 0.0D0
    call MPI_SENDRECV(cden(nxs1-1),2,mnpr,ndown,200, &
                      ctmp(1),2,mnpr,nup,200,        &
                      ncomw,nstat,nerr)
    cden(nxe-1) = cden(nxe-1)+ctmp(1)
    cden(nxe)   = cden(nxe)  +ctmp(2)

    ctmp(1:2) = 0.0D0
    call MPI_SENDRECV(cden(nxe1),2,mnpr,nup,100, &
                      ctmp(1),2,mnpr,ndown,100,    &
                      ncomw,nstat,nerr)
    cden(nxs)   = cden(nxs)  +ctmp(1)
    cden(nxs+1) = cden(nxs+1)+ctmp(2)

    if(bc == -1)then
       !reflective condition
       if(nrank == nroot)then
          cden(nxs)   = cden(nxs)  -cden(nxs1)
          cden(nxs+1) = cden(nxs+1)-cden(nxs1-1)
       endif
       if(nrank == nproc-1)then
          cden(nxe-2) = cden(nxe-2)-cden(nxe1)
          cden(nxe-1) = cden(nxe-1)-cden(nxe)
       endif
    endif

  end subroutine boundary__charge



end module boundary
