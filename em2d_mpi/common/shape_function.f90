module shape_function

  implicit none

  private

  public :: sf

contains 

  function sf(ip,rp,nsfo)

    integer, intent(in) :: ip, nsfo
    real(8), intent(in) :: rp
    real(8) :: sf(5)
    real(8) :: dx

    select case(nsfo)

    !NGP
    case(0)
       sf(1) = 0.D0
       sf(2) = 0.D0
       sf(3) = 1.D0
       sf(4) = 0.D0
       sf(5) = 0.D0

    !CIC
    case(1)
       dx = rp-ip

       sf(1) = 0.D0
       sf(2) = 0.D0
       sf(3) = 1.-abs(dx)
       sf(4) = 0.D0
       sf(5) = 0.D0

       sf(3+int(sign(1.D0,dx))) = abs(dx)

    !spline
    case(2)
       dx = rp-ip
       sf(1) = 0.D0
       sf(2) = 0.5*(0.5-dx)**2
       sf(3) = 0.75-dx*dx
       sf(4) = 0.5*(0.5+dx)**2
       sf(5) = 0.D0

    !CIC is the default shape function
    case default
       dx = rp-ip

       sf(1) = 0.D0
       sf(2) = 0.D0
       sf(3) = 1.-abs(dx)
       sf(4) = 0.D0
       sf(5) = 0.D0

       sf(3+int(sign(1.D0,dx))) = abs(dx)

    end select

  end function sf

end module shape_function
