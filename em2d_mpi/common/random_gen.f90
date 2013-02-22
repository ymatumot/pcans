module random_gen

  implicit none

  private

  public :: random_gen__init, random_gen__bm, random_gen__sobol, random_gen__LLT

  real(8), save :: pi


contains


  subroutine random_gen__init(nrank)

    integer, intent(in)  :: nrank
    integer              :: n
    integer, allocatable :: seed(:)

    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)+nrank
    call random_seed(put=seed)
    deallocate(seed)

    pi = 4.*atan(1.D0)

  end subroutine random_gen__init


  subroutine random_gen__bm(r1,r2)
    !Box-Muller transform

    real(8), intent(out) :: r1, r2
    real(8)              :: aa, bb

    call random_number(aa) ! [0,1)
    call random_number(bb)
    aa = 1-aa ! (0,1]

    r1 = dsqrt(-2.*dlog(aa))*cos(2.*pi*bb)
    r2 = dsqrt(-2.*dlog(aa))*sin(2.*pi*bb)

  end subroutine random_gen__bm


  subroutine random_gen__sobol(u1,u2,u3,T)
    !! Sobol's algorithm

    real(8), intent(in)  :: T
    real(8), intent(out) :: u1, u2, u3
    real(8)              :: r1, r2, r3, r4
    real(8)              :: e3, e4
    real(8)              :: aa, bb, cc

    do
       !! Gamma distributions
       call random_number(r1) ! [0,1)
       call random_number(r2) ! [0,1)
       call random_number(r3) ! [0,1)
       call random_number(r4) ! [0,1)
       r1 = 1-r1 ! (0,1]
       r2 = 1-r2 ! (0,1]
       r3 = 1-r3 ! (0,1]
       r4 = 1-r4 ! (0,1]
       e3 = - T * dlog(r1*r2*r3)
       e4 = - T * dlog(r1*r2*r3*r4)
       !! criterion
       if( ( e4**2 - e3**2 ) .gt. 1.d0 ) then
          exit
       endif
    enddo

    !! Spherical scattering
    call random_number(aa)
    call random_number(bb)
    call random_number(cc)

    aa = 2.d0 * aa - 1.d0
    bb = dsqrt( 1.d0 - aa**2 )
    u1 = e3 * aa
    u2 = e3 * bb * cos(2.*pi*cc)
    u3 = e3 * bb * sin(2.*pi*cc)

  end subroutine random_gen__sobol


  subroutine random_gen__LLT(u1,u2,u3,ubulk)
    !! Lossless Lorentz Transformation (specific case only)

    real(8), intent(in)    :: u1, u2, ubulk
    real(8), intent(inout) :: u3
    real(8)                :: u0, dice
    real(8)                :: beta, bgamma

    bgamma = dsqrt( 1.d0 + ubulk**2 )
    beta   = ubulk / bgamma
    u0     = dsqrt(1.d0+u1**2+u2**2+u3**2)
    if( (u3*beta) .lt. 0.d0 ) then
       call random_number(dice)
       if( (u0*dice).lt.( -beta*u3 ) )then
          u3 = - u3
       endif
    endif
    u3 = bgamma*u3 + u0*ubulk
       
  end subroutine random_gen__LLT


end module random_gen

