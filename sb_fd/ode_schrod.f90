MODULE ode_schrod
  IMPLICIT NONE
  !=============================================================================!
  ! This module contains:                                                       !
  ! read_pot   : the subroutine to read the potential and energy from the file  !
  ! f(x)       : a function to define the ODE (schrodinger eq)                  !
  ! v(x)       : a function to define the potential at each x                   !
  ! find_coeff : a subroutine to calculate the T and R coefficients             !
  ! gauss()    : a function to find the gaussian distribution of energies       !
  !=============================================================================!


  !===============================DECLARATIONS==================================!
  !constants
  integer, parameter             :: dp=selected_real_kind(15,300)
  real(kind=dp), parameter       :: pi        = 3.141592653589793238462643_dp
  !dimensions and shape of the potential
  character(len=20), save        :: pot_form
  real(kind=dp), save            :: a,b,c,d,height,E, e_start, e_finish, Tr, refl
  integer, save                  :: e_it, num_e
  real(kind=dp), save            :: para_cons, temperature
  !variables
  integer                        :: ierr


contains

  subroutine read_pot
    IMPLICIT NONE
    !===========================================================================!
    ! Reads the form and dimentions of the given potential from a file.         !
    ! Then asks the user for values of e_start, e_finish and e_it               !
    !---------------------------------------------------------------------------!
    ! Input:                                                                    !
    ! input2.dat                                                                !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! pot_form, height, a, b, c, d, para_cons, e_start, e_finish, e_it          !
    !===========================================================================!

    !Open the input file to read the data from
    open(10,file='input2.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening input.dat"

    !Read the potential data
    read(10,*) pot_form
    if (pot_form=='triangle') then
       read(10,*) height,a,b
    else if (pot_form=='parabolic_triangle') then
       read(10,*) height,a,b
    else if (pot_form=='rectangle') then
       read(10,*) height,a,b
    else if (pot_form=='double_rectangle') then
       read(10,*) height,a,b,c,d
    else if (pot_form=='trahttps://github.com/nf679/sb_et.gitpezoid') then
       read(10,*) height,a,b,c,d
    else
       stop "Invalid potential name"
    end if

    !This is calculated for the parabolic triangle case (ignore otherwise)
    para_cons = height/(a**2 - 2.0_dp*b*a + b**2)

    !Read the energy values from the file
    read(10,*) e_start
    read(10,*) e_finish
    read(10,*) e_it

    !Read the number of particles
    read(10,*) num_e

    !Read the temperature to calculate the Gaussian
    read(10,*) temperature
    !Put the temperature in atomic units (make this a function at some point)
    temperature = temperature / 3.1577464e5_dp

    !Set the value of E to the starting energy
    E = e_start

    !Close the input file
    close(10,iostat=ierr)
    if(ierr/=0) stop "Error in closing input.dat"

    return
  end subroutine read_pot

  function f(x,w)
    IMPLICIT NONE
    !===========================================================================!
    ! Defines the ode under study. This is the Schrödinger Eqn.                 !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! function V(x)                                                             !
    !---------------------------------------------------------------------------!
    ! Input:                                                                    !
    ! x  : Distance (Bohrs)                                                     !
    ! w  : Matrix containing psi and psi' (Hartrees)                            !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! f (matrix containing the ode)                                             !
    !===========================================================================!


    !===============================DECLARATIONS==================================!
    real(kind=dp), intent(in)                  :: x
    complex(kind=dp), dimension(2), intent(in) :: w
    complex(kind=dp), dimension(2)             :: f

    !Set the components of the ODE
    f(1) = w(2)
    f(2) = 2*(V(x)-E)*w(1)

    return
  end function f

  real(kind=dp) function V(x)
    IMPLICIT NONE
    !===========================================================================!
    ! This function take the constants read in by read_pot and THEN             !
    ! calculates a value for the potential at a particular value of x           !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! read_pot                                                                  !
    !---------------------------------------------------------------------------!
    ! Input:                                                                    !
    ! x (Bohrs)                                                                 !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! V (single value of potential)                                             !
    !===========================================================================!


    !===============================DECLARATIONS==================================!
    real(kind=dp), intent(in) :: x

    !Triangle potential
    if (pot_form=='triangle') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V= height*x/(a-b) + height*b/(b-a)
       else if (x .gt. b) then
          V=0.0_dp
       end if

    !Parabolic traingle potential
    else if (pot_form=='parabolic_triangle') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V= para_cons*x**2 - 2.0_dp*b*para_cons*x + para_cons*b**2
       else if (x .gt. b) then
          V=0.0_dp
       end if

    !Rectangle potential
    else if (pot_form=='rectangle') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V=height
       else if (x .gt. b) then
          V=0.0_dp
       end if

    !Double rectangle potential
    else if (pot_form=='double_rectangle') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V=height
       else if (x .gt. b .and. x .lt. c) then
          V=0.0_dp
       else if (x .ge. c .and. x .le. d) then
          V=height
       else if (x .gt. d) then
          V=0.0_dp
       end if

    !Trapezoid potential
    else if (pot_form=='trapezoid') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V= height*(x-a)/(b-a)
       else if (x .gt. b .and. x .lt. c) then
          V=height
       else if (x .ge. c .and. x .le. d) then
          V= height + height*(x-c)/(c-d)
       else if (x .gt. d) then
          V=0.0_dp
       end if

    !For any other shape of potential
    else
       stop "Error finding potential for the input shape"
    end if

    return
  end function V


  function gauss(E)
    IMPLICIT NONE
    !===========================================================================!
    ! A function to calculate the gaussian spread of energies of the system     !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! height, temperature, pi                                              !
    !===========================================================================!

    !=======================DECLARE VARIABLES====================================!
    real(kind=dp) :: gauss
    real(kind=dp), intent(in) :: E

    !Equation for the gaussian distribution
    gauss = (1.0_dp / (sqrt(2.0_dp*pi*(temperature**2)))) * exp(((-1.0_dp*(E - height)**2) / (2.0_dp * temperature**2)))

  end function gauss

  function fermi_d(E)
    IMPLICIT NONE
    !============================================================================!
    ! A function to calculate the fermi dirac spread of energies of the system   !
    !----------------------------------------------------------------------------!
    ! Requires:                                                                  !
    ! height, temperature                                                        !
    !============================================================================!

    !=========================DECLARE VARIABLES==================================!
    real(kind=dp) :: fermi_d
    real(kind=dp), intent(in) :: E

    !Equation for the fermi dirac distribution
    fermi_d = 1.0_dp/((exp((E-height)/temperature)) + 1.0_dp)


  subroutine find_coeff
    IMPLICIT NONE
    !===========================================================================!
    ! This subroutine goes through the newly created output.dat file and gets   !
    ! the values of psi at x=0 and x=pi/2*k. Then, the values are used to       !
    ! calculate the reflection and transmission coefficients.                   !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! output.dat         : data file created by runge_kutta subroutine          !
    !---------------------------------------------------------------------------!                                                                          !
    ! Input:                                                                    !
    !  E                 : Energy (Hartrees)                                    !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! T, R               : Transmission / Reflection coefficients               !
    ! modeltest.dat      : data file containing the P+S E vs T data             !
    ! Energy_v_coeff.dat : data file containing calculated E vs T data          !
    !===========================================================================!


    !===============================DECLARATIONS==================================!
    complex(kind=dp) :: psi_1,psi_2,At,Ar
    real(kind=dp)    :: x,Re,Im,k, sc_parameter, delta_E
    integer          :: ierr

    !Set the value of the parameter for the padovanni and stratton model
    sc_parameter = 1.0e-5_dp

    !Open the file to read the wavefunction from
    open(15,file='output.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening output.dat"

    !Find psi at x = 0
    read(15,*) x,Re,Im
    loop1 : do
       if (abs(x) .gt. 1.0e-6_dp) then
          read(15,*) x,Re,Im
       else
          exit loop1
       end if
    end do loop1

    !Store the value
    psi_1 = cmplx(Re,Im,kind=dp)

    !Close the file
    close(15,iostat=ierr)
    if (ierr/=0) stop "Error in close output.dat"

    !Open the file again to read the next wavefunction from
    open(15,file='output.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening output.dat"

    !Find psi at x = pi/2*k
    k = sqrt(2*E)
    read(15,*) x,Re,Im
    loop2 : do
       if (abs(x) .gt. pi*0.5_dp/k) then
          read(15,*) x,Re,Im
       else
          exit loop2
       end if
    end do loop2

    !Store the value
    psi_2 = cmplx(Re,Im,kind=dp)

    !Close the file
    close(15,iostat=ierr)
    if (ierr/=0) stop "Error in close output.dat"

    !Calculate Reflection coeff
    At   = (psi_1 - cmplx(0.0_dp,1.0_dp,kind=dp)*psi_2)*0.5_dp
    Ar   = (cmplx(0.0_dp,1.0_dp,kind=dp)*psi_2 + psi_1)*0.5_dp
    refl = (abs(Ar)**2)/(abs(At)**2)

    !Test how this changes the curve
    !if(E >= height) then
    !refl = refl*(-1.0_dp)
    !end if

    !Calculate T using Padovanni and Stratton Formula
    delta_E = height - E
    if (delta_E < 0.0_dp) then
      delta_E = delta_E*(-1.0_dp)
    end if
    Tr = exp((-(2.0_dp/3.0_dp)*(delta_E)**(3.0_dp/2.0_dp))/(sc_parameter))


    return
  end subroutine find_coeff

end module ode_schrod
