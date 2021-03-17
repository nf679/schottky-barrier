module rk
  
  use ode_schrod
  implicit none
  
contains
  
  subroutine runge_kutta
    !===========================================================================!
    ! This is an algorithm for the fourth order Runge Kutta method. It is a 4th !
    ! order ODE solver. However, this method is written so that it starts on    !
    ! b and moves backwards, i.e., opposite to the normal direction             !
    !---------------------------------------------------------------------------!
    ! Input:                                                                    !
    !                                                                           !
    !  a,b : start and end points of the system                                 !
    !  n   : number of iterations (also defines step size through h)            !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    !  ode_schrod(or test) must be called before running this subroutine        !
    !---------------------------------------------------------------------------!
    ! This input is taken from a file called input.dat, which must be present   !
    ! when the program is run                                                   !
    !===========================================================================!
    implicit none
    
    !parameters
    integer                        :: i,n,ierr
    real(kind=dp)                  :: x,h,start,finish,wave_n
    complex(kind=dp), dimension(2) :: w,k1,k2,k3,k4
    
    !read in input parameters
    open(10,file='input1.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening input1.dat"
    read(10,*) start,finish
    read(10,*) n
    close(unit=10,iostat=ierr)
    if (ierr/=0) stop "Error in closing input1.dat"

    !calculate starting psi and psi'
    wave_n = sqrt(2.0_dp*E)
    w(1)  = exp(cmplx(0.0_dp,wave_n,kind=dp)*finish)
    w(2)  = cmplx(0.0_dp,wave_n,kind=dp)*w(1)        
    h     = (finish-start)/real(n,dp)
    x     = finish


    !Start writing output data
    open(15,file='output.dat',status='replace',iostat=ierr)
    if (ierr/=0) stop "Error in opening output.dat"
    write(15,*) x,real(w(1)),aimag(w(1)),real(w(2)),aimag(w(2))
    
    do i=1,n
       ! 4 function calls, each depends on the preceding one
       k1 = h*f(x,w)
       k2 = h*f(x-(h/2.0_dp),w-(k1/2.0_dp))
       k3 = h*f(x-(h/2.0_dp),w-(k2/2.0_dp))       
       k4 = h*f(x-h,w-k3)
       
       w  = w - (k1+2.0_dp*k2+2.0_dp*k3+k4)/6.0_dp
       x  = x - h
       write(15,*) x,real(w(1)),aimag(w(1)),real(w(2)),aimag(w(2))
       
    end do

    close(15,iostat=ierr)
    if (ierr/=0) stop "Error in closing output.dat"
    
    return
  end subroutine runge_kutta
        
end module rk


program ode_solver
  use rk  
  implicit none
  integer :: i, ierr
  real(kind=dp) :: e_h

  call read_pot
  e_h = (e_finish - e_start) / real(e_it)

  do i=0,(e_it-1)
    call runge_kutta
    call find_coeff

    E = E + e_h
  end do
  

end program ode_solver
