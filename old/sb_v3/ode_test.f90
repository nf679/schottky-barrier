module ode_test
  implicit none
  
  integer, parameter             :: dp=selected_real_kind(15,300)
  real(kind=dp),save             :: E = 2.0_dp
  
contains
  
  function f(x,w)
    !========================================================================!
    ! This suroutine is an analytically solvable ode that can be used as a   !
    ! test. (It is the most basic case of the simple harmonic oscillator)    !
    !========================================================================!
    implicit none
    
    complex(kind=dp), dimension(2)             :: f
    complex(kind=dp), dimension(2), intent(in) :: w
    real(kind=dp)                              :: x 
    
    f(1) =  w(2)
    f(2) = -w(1)    
    
    return
  end function f
  
end module ode_test
