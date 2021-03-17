module ode_schrod
  implicit none
  
  integer, parameter             :: dp=selected_real_kind(15,300)

  !constants
  real(kind=dp), parameter :: pi        = 3.141592653589793238462643_dp

  !dimensions and shape of the potential
  character(len=20), save        :: pot_form
  real(kind=dp), save            :: a,b,c,d,height,E, e_start, e_finish
  integer, save                  :: e_it
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

    close(10,iostat=ierr)
    if(ierr/=0) stop "Error in closing input.dat"
  

    print*, "What do you want as your starting energy?(Real number)"
    read(*,*) e_start
    E = e_start
    print*, "What do you want your final energy to be? (Real number)"
    read(*,*) e_finish
    print*, "Number of iterations? (Integer)"
    read(*,*) e_it

    
    return
  end subroutine read_pot
  
  function f(x,w)
    !======================================================================!
    ! Defines the ode under study. This is the Schrödinger Eqn             !
    !----------------------------------------------------------------------!
    ! Requires function V(x)                                               !
    !----------------------------------------------------------------------!
    ! Units:                                                               !
    ! The TISE has been written using atomic units, therefore              !
    ! x is in Hartrees                                                     !
    ! E is in Bohrs                                                        !
    !======================================================================!
    implicit none
    real(kind=dp), intent(in)                  :: x
    complex(kind=dp), dimension(2), intent(in) :: w
    complex(kind=dp), dimension(2)             :: f
    
    f(1) = w(2)
    f(2) = 2*(V(x)-E)*w(1)
    
    return
  end function f
  
  real(kind=dp) function V(x)
    !======================================================================!
    ! This function take the constants read in by read_pot and then        !
    ! calculates a value for the potential at a particular value of x      !
    !----------------------------------------------------------------------!
    ! Requires:                                                            !
    ! read_pot                                                             !
    !----------------------------------------------------------------------!
    ! x is in units of Bohr (atomic units)                                 !
    !======================================================================!
    implicit none
    real(kind=dp), intent(in) :: x
    
    if (pot_form=='triangle') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V= height*x/(a-b) + height*b/(b-a)
       else if (x .gt. b) then
          V=0.0_dp
       end if

    else if (pot_form=='parabolic_triangle') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V= para_cons*x**2 - 2.0_dp*b*para_cons*x + para_cons*b**2
       else if (x .gt. b) then
          V=0.0_dp
       end if
       
    else if (pot_form=='rectangle') then
       if (x .lt. a) then
          V=0.0_dp
       else if (x .ge. a .and. x .le. b) then
          V=height
       else if (x .gt. b) then
          V=0.0_dp
       end if
       
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
    else
       stop "Error finding potential for the input shape"
    end if
    
    return
  end function V

  subroutine find_coeff
    !==========================================================================!
    ! This subroutine goes through the newly created output.dat file and gets  !
    ! the values of psi at x=0 and x=pi/2*k. Then, the values are used to      !
    ! calculate the Reflection coefficient (and hence Transmission)            !
    !--------------------------------------------------------------------------!
    ! Subroutine must occur after the runge_kutta subroutine.                  !
    ! This is becuase it required the file output.dat to exist                 !
    !--------------------------------------------------------------------------!
    !                                                                          !
    ! Input:                                                                   !
    !  E : Energy (Hartrees)                                                   !     
    !==========================================================================!
    implicit none
    
    complex(kind=dp) :: psi_1,psi_2,At,Ar
    real(kind=dp)    :: x,Re,Im,k,refl, Tr, sc_parameter, delta_E
    integer          :: ierr
    
    sc_parameter = 1.0e-5_dp

    open(15,file='output.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening output.dat"
    
    !find psi at x = 0
    read(15,*) x,Re,Im
    loop1 : do
       if (abs(x) .gt. 1.0e-6_dp) then
          read(15,*) x,Re,Im
       else
          exit loop1
       end if
    end do loop1

    !store value
    psi_1 = cmplx(Re,Im,kind=dp)
    
    close(15,iostat=ierr)
    if (ierr/=0) stop "Error in close output.dat"
    open(15,file='output.dat',status='old',iostat=ierr)
    if (ierr/=0) stop "Error in opening output.dat"
    
    !find psi at x = pi/2*k
    k = sqrt(2*E)
    read(15,*) x,Re,Im
    loop2 : do
       if (abs(x) .gt. pi*0.5_dp/k) then
          read(15,*) x,Re,Im
       else
          exit loop2
       end if
    end do loop2
    
    !store value
    psi_2 = cmplx(Re,Im,kind=dp)
    
    close(15,iostat=ierr)
    if (ierr/=0) stop "Error in close output.dat"

    !calculate Reflectoin coeff
    At   = (psi_1 - cmplx(0.0_dp,1.0_dp,kind=dp)*psi_2)*0.5_dp
    Ar   = (cmplx(0.0_dp,1.0_dp,kind=dp)*psi_2 + psi_1)*0.5_dp    
    refl = (abs(Ar)**2)/(abs(At)**2)

   ! if(E >= height) then 
     ! refl = refl*(-1.0_dp)
  !  end if

    !calculate t using p+s
    delta_E = height - E
    if (delta_E < 0.0_dp) then
      delta_E = delta_E*(-1.0_dp)
    end if
    Tr = exp((-(2.0_dp/3.0_dp)*(delta_E)**(3.0_dp/2.0_dp))/(sc_parameter))

    open(30,file='modeltest.dat',position='append',iostat=ierr)
    write(30,*) E, Tr
    close(30,iostat=ierr)
    


    !!output energy and coefficients
    open(25,file='Energy_v_Coeff.dat',position='append',iostat=ierr)
    if (ierr/=0) stop "Error in opening Energy_v_Coeff.dat"
    !

      write(25,*) E,(1.0_dp-refl)


    !
    close(20,iostat=ierr)
    if (ierr/=0) stop "Error closing Energy_v_Coeff.dat"
    

    return
  end subroutine find_coeff
    
end module ode_schrod
