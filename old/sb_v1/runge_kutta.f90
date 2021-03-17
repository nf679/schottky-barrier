module rk
  use ode_schrod
  implicit none
  !=============================================================================!
  ! This module contains:                                                       !
  ! runge_kutta         : subroutine to calculate psi and psi'                  !
  !=============================================================================!
  
contains
  
  subroutine runge_kutta
    implicit none 
    !===========================================================================!
    ! This is an algorithm for the fourth order Runge Kutta method.             !
    ! This method is written so that it starts on b and moves backwards.        !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! ode_schrod must be called before running this subroutine                  !
    !---------------------------------------------------------------------------!
    ! Input:                                                                    !
    ! input1.dat : data file containing start, finish, n                        !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! output1.dat : data file containing x, psi                                 !
    ! OR                                                                        !
    ! output2.dat : data file containing x, psi'                                !
    !===========================================================================!
    
    
    !===============================DECLARATIONS==================================!
    integer                        :: i,n,ierr
    real(kind=dp)                  :: x,h,start,finish,wave_n
    complex(kind=dp), dimension(2) :: w,k1,k2,k3,k4
    
    !Open the file to read the input parameters from
    open(10,file='input1.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening input1.dat"

    !Read the input parameters
    read(10,*) start,finish
    read(10,*) n

    !Close the input file
    close(unit=10,iostat=ierr)
    if (ierr/=0) stop "Error in closing input1.dat"

    !Calculate starting psi and psi'
    wave_n = sqrt(2.0_dp*E)
    w(1)  = exp(cmplx(0.0_dp,wave_n,kind=dp)*finish)
    w(2)  = cmplx(0.0_dp,wave_n,kind=dp)*w(1)        
    h     = (finish-start)/real(n,dp)
    x     = finish

    !Start writing output data for option 1, 3 and 4
    if (output_option == 1 .or. output_option == 3 .or. output_option == 4) then
      open(15,file='output1.dat',status='replace',iostat=ierr)
      if (ierr/=0) stop "Error in opening output1.dat"
      write(15,*) x,real(w(1)),aimag(w(1))

    !Start writing output data for option 2
    else if (output_option == 2) then
      open(16,file='output2.dat',status='replace',iostat=ierr)
      if (ierr/=0) stop "Error in opening output2.dat"
      write(16,*) x,real(w(2)),aimag(w(2))
    end if


    !Do loop for the RK4 method
    do i=1,n
       k1 = h*f(x,w)
       k2 = h*f(x-(h/2.0_dp),w-(k1/2.0_dp))
       k3 = h*f(x-(h/2.0_dp),w-(k2/2.0_dp))       
       k4 = h*f(x-h,w-k3)
       w  = w - (k1+2.0_dp*k2+2.0_dp*k3+k4)/6.0_dp
       x  = x - h

       !For option 1,3,4 write the data for psi to the file
       if (output_option == 1 .or. output_option == 3 .or. output_option == 4) then
         write(15,*) x,real(w(1)),aimag(w(1))

       !For option 2 write the data for psi' to the file   
       else if (output_option == 2) then
         write(16,*) x,real(w(2)),aimag(w(2))
       end if
       
    end do

    !Close output1.dat
    if (output_option == 1 .or. output_option == 3 .or. output_option == 4) then
      close(15,iostat=ierr)
      if (ierr/=0) stop "Error in closing output1.dat"
    end if

    !Close output2.dat
    if (output_option == 2) then
      close(16,iostat=ierr)
      if (ierr/=0) stop "Error in closing output2.dat"
    end if
    
    return
  end subroutine runge_kutta
        
end module rk


program ode_solver
  use rk  
  implicit none
  !=============================================================================!
  ! The main program which calls all of the subroutines and does different      !
  ! things depending on the option selected (See README file).                  !
  !=============================================================================!


  !===============================DECLARATIONS==================================!
  integer :: i
  real(kind=dp) :: E_h
  
  !Call the subroutines to read the potential data and calculate psi using RK4
  call read_pot
  call runge_kutta

  !For option 3 and 4 write data for E and T to a file
  if (output_option == 3 .or. output_option == 4) then
    !Calculate the step size for the energy calculation
    E_h = ((E_start) - (E_finish)) / real(E_n,dp)
    !Loop over the specified energy range
    do i=0, E_n
      E = E_start + real(i*E_h)
      call runge_kutta
      call find_coeff
    end do
  end if

  
end program ode_solver
