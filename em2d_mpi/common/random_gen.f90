module random_gen

  implicit none

  private

  public :: random_gen__init, random_gen__bm, random_gen__sobol

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

    call random_number(aa)
    call random_number(bb)

    r1 = dsqrt(-2.*dlog(aa))*cos(2.*pi*bb)
    r2 = dsqrt(-2.*dlog(aa))*sin(2.*pi*bb)

  end subroutine random_gen__bm


  subroutine random_gen__sobol(T,u1,u2,u3)
    !! Sobol's algorithm

    real(8), intent(in)  :: T
    real(8), intent(out) :: u1, u2, u3
    real(8)              :: r1, r2, r3, r4
    real(8)              :: e3, e4
    real(8)              :: aa, bb, cc

    do
       !! Gamma distributions
       call random_number(r1)
       call random_number(r2)
       call random_number(r3)
       call random_number(r4)
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


end module random_gen

