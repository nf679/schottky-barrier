module ode_schrod
  implicit none
  
  integer, parameter             :: dp=selected_real_kind(15,300)

  !constants
  real(kind=dp), parameter :: pi        = 3.141592653589793238462643_dp

  !dimensions and shape of the potential
  character(len=20), save        :: pot_form
  integer, save                  :: output_option
  real(kind=dp), save            :: a,b,c,d,height
  real(kind=dp),save             :: E, E_start, E_finish 
  integer                        :: E_n 
  real(kind=dp), save            :: para_cons
  
  
contains
  
  subroutine read_pot
    !======================================================================!
    ! Reads the form and dimentions of the given potential from a file,    !
    ! input.dat (This will already have been opened previously in a        !
    ! different subroutine)                                                !
    !----------------------------------------------------------------------!
    ! This subroutine need to be called before using the Runge Kutta       !
    ! method, if the ode is giong to be the shrodinger equation            !
    !======================================================================!
    implicit none
    
    integer :: ierr
    
    open(10,file='input2.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening input.dat"

    ! read the output type and make sure its one of the correct options
    read(10,*) output_option
    if (output_option .ne. 1 .and. output_option .ne. 2 .and. output_option &
    .ne. 3 .and. output_option .ne. 4 .and. output_option .ne. 5) then
      stop "Invalid output option, please enter a number between 1 and 5"
    end if 

    read(10,*) pot_form
    if (pot_form=='triangle') then
       read(10,*) height,a,b
    else if (pot_form=='parabolic_triangle') then
       read(10,*) height,a,b
    else if (pot_form=='rectangle') then
       read(10,*) height,a,b
    else if (pot_form=='double_rectangle') then
       read(10,*) height,a,b,c,d
    else if (pot_form=='trapezoid') then
       read(10,*) height,a,b,c,d
    else
       stop "Invalid potential name"
    end if

    !this is calculated for the parabolic triangle case (ignore otherwise)
    para_cons = height/(a**2 - 2.0_dp*b*a + b**2)
    
    read(10,*) E

    !if option 3 or 4 is selected, we need the start and finish energies
    if (output_option == 3 .or. output_option == 4) then
      read(10,*) E_start
      read(10,*) E_finish
      read(10,*) E_n
    end if

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
  
  function V(x)
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
    real(kind=dp) :: V
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


  subroutine find_coeff
    IMPLICIT NONE
    !===========================================================================!
    ! This subroutine goes through the newly created output.dat file and gets   !
    ! the values of psi at x=0 and x=pi/2*k. Then, the values are used to       !
    ! calculate the reflection and transmission coefficients.                   !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! output1.dat        : data file created by runge_kutta subroutine          !
    !---------------------------------------------------------------------------!                                                                          !
    ! Input:                                                                    !
    ! E                  : Energy (Hartrees)                                    ! 
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! T, R               : Transmission / Reflection coefficients               !  
    ! outpu3.dat : data file containing calculated E vs T data                  !
    ! OR                                                                        !
    ! output4.dat : data file containing calculated E vs R data                 !
    !===========================================================================!


    !===============================DECLARATIONS==================================!
    complex(kind=dp) :: psi_1,psi_2,At,Ar
    real(kind=dp)    :: x,Re,Im,k,refl
    integer          :: ierr
    
    !Open a file to read the wavefunction data from
    open(15,file='output1.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening output1.dat"
    
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
    if (ierr/=0) stop "Error in close output1.dat"

    !Reopen the file for the next value of psi we want
    open(15,file='output1.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening output1.dat"
    
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
    if (ierr/=0) stop "Error in close output1.dat"

    !Calculate Reflection coeff
    At   = (psi_1 - cmplx(0.0_dp,1.0_dp,kind=dp)*psi_2)*0.5_dp
    Ar   = (cmplx(0.0_dp,1.0_dp,kind=dp)*psi_2 + psi_1)*0.5_dp    
    refl = (abs(Ar)**2)/(abs(At)**2)

    !For option 3 write the E,T data to a file
    if (output_option == 3) then        
      open(20,file='output3.dat',position='append',iostat=ierr)
      if (ierr/=0) stop "Error in opening output3.dat"
      write(20,*) E,1.0_dp-refl
      close(20,iostat=ierr)
      if (ierr/=0) stop "Error closing output3.dat"

    !For option 4, write the E,R data to a file
    else if (output_option == 4) then
      open(25,file='output4.dat',position='append',iostat=ierr)
      if (ierr/=0) stop "Error in opening output4.dat"
      write(25,*) E,refl
      close(25,iostat=ierr)
      if (ierr/=0) stop "Error closing output4.dat"
    end if
    
    !For the other options, print the value of T to the console
    if (output_option .ne. 3 .and. output_option .ne. 4) then
      print*, 1.0_dp-refl
    end if

    return
  end subroutine find_coeff
    
end module ode_schrod
