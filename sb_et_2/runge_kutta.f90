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
    ! output.dat : data file containing x, psi, psi'                            !
    !===========================================================================!

    !===============================DECLARATIONS==================================!
    integer                        :: i,n,ierr
    real(kind=dp)                  :: x,h,start,finish,wave_n
    complex(kind=dp), dimension(2) :: w,k1,k2,k3,k4

    !Open a file to read the input parameters
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

    !Start writing output data
    open(15,file='output.dat',status='replace',iostat=ierr)
    if (ierr/=0) stop "Error in opening output.dat"
    write(15,*) x,real(w(1)),aimag(w(1)),real(w(2)),aimag(w(2))

    !Do loop to calculate psi and psi' at each x using RK4
    do i=1,n
       k1 = h*f(x,w)
       k2 = h*f(x-(h/2.0_dp),w-(k1/2.0_dp))
       k3 = h*f(x-(h/2.0_dp),w-(k2/2.0_dp))
       k4 = h*f(x-h,w-k3)
       w  = w - (k1+2.0_dp*k2+2.0_dp*k3+k4)/6.0_dp
       x  = x - h

      !Write the data to the file
       write(15,*) x,real(w(1)),aimag(w(1)),real(w(2)),aimag(w(2))

    end do

    !Close the output file
    close(15,iostat=ierr)
    if (ierr/=0) stop "Error in closing output.dat"

    return
  end subroutine runge_kutta

end module rk


program ode_solver
  use rk
  implicit none
  !=============================================================================!
  ! The main program which calls all of the subroutines, gets user input of     !
  ! energy range to loop over, and number of iterations. It then loops over     !
  ! this energy range to calculate T for each energy value.                     !
  !=============================================================================!

  !===============================DECLARATIONS==================================!
  integer :: i
  real(kind=dp) :: e_h

  !Read the potential data
  call read_pot

  !Calculate the step size for the energy loop
  e_h = (e_finish - e_start) / real(e_it)

  !Open the files to write the data to
  open(30,file='model.dat',status='replace',iostat=ierr)
  if (ierr/=0) stop "Error in opening model.dat"
  open(25,file='Energy.dat',status='replace',iostat=ierr)
  if (ierr/=0) stop "Error in opening Energy.dat"


  !Loop over the energy range to get E v T data
  do i=0,(e_it-1)

    call runge_kutta
    call find_coeff

    !Write the P+S data to the file
    write(30,*) E, Tr

    !Write the calculated E vs T data to the file
    write(25,*) E,(1.0_dp-refl)

    !Increment the energy 
    E = E + e_h

  end do

  !Close the files
  close(30,iostat=ierr)
  if (ierr/=0) stop "Error closing model.dat"
  close(25,iostat=ierr)
  if (ierr/=0) stop "Error closing Energy.dat"

end program ode_solver
