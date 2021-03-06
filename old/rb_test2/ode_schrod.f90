MODULE ode_schrod
  IMPLICIT NONE
  !=============================================================================!
  ! This module contains:                                                       !
  ! read_pot   : the subroutine to read the potential and energy from the file  !
  ! f(x)       : a function to define the ODE (schrodinger eq)                  !
  ! v(x)       : a function to define the potential at each x                   !
  ! find_coeff : a subroutine to calculate the T and R coefficients             !
  !=============================================================================!


  !===============================DECLARATIONS==================================!
  !Constants
  INTEGER, PARAMETER             :: dp=selected_real_kind(15,300)
  REAL(kind=dp), PARAMETER       :: pi        = 3.141592653589793238462643_dp
  !Dimensions and shape of the potential
  CHARACTER(len=20), SAVE        :: pot_form
  REAL(kind=dp), SAVE            :: a,b,c,d,height,E, bias_i, bias_f, bias_h
  REAL(kind=dp), SAVE            :: bias 
  REAL(kind=dp), SAVE            :: para_cons, refl 
  INTEGER, SAVE                  :: bias_n
  !Variables
  INTEGER                        :: ierr

CONTAINS
  
  SUBROUTINE read_pot
    IMPLICIT NONE
    !===========================================================================!
    ! Reads the form and dimentions of the given potential from a file.         !
    ! THEN reads the single value of E and bias.                                !
    !---------------------------------------------------------------------------!
    ! Input:                                                                    !
    ! input2.dat                                                                !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! pot_form, height, a, b, c, d, para_cons, E, bias                          !
    !===========================================================================!
    
    !Open the file
    OPEN(20,file='input2.dat',status='old',iostat=ierr)
    IF (ierr/=0) STOP "Error in opening input2.dat"

    !Read the potential parameters
    READ(20,*) pot_form
    IF (pot_form=='triangle') THEN
       READ(20,*) height,a,b
    ELSE IF (pot_form=='parabolic_triangle') THEN
       READ(20,*) height,a,b
    ELSE IF (pot_form=='rectangle') THEN
       READ(20,*) height,a,b
    ELSE IF (pot_form=='double_rectangle') THEN
       READ(20,*) height,a,b,c,d
    ELSE IF (pot_form=='trapezoid') THEN
       READ(20,*) height,a,b,c,d
    ELSE
       STOP "Invalid potential name"
    END IF
    
    !Read the value of E from the file
    READ(20,*) E

    !Read the bias parameters from the file
    READ(20,*) bias_i
    READ(20,*) bias_f
    READ(20,*) bias_n

    !Calculate the step size for the bias parameters
    bias_h = (bias_f - bias_i) / REAL(bias_n)

    !Set the initial value of bias to bias_i#
    bias = bias_i

    !Close the file
    CLOSE(20,iostat=ierr)
    IF(ierr/=0) STOP "Error in closing input2.dat"
    
    RETURN
  END SUBROUTINE read_pot

  
  FUNCTION f(x,w)
    IMPLICIT NONE
    !===========================================================================!
    ! Defines the ode under study. This is the Schr??dinger Eqn.                 !
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


    !==============================DECLARATIONS=================================!
    REAL(kind=dp), INTENT(in)                  :: x
    COMPLEX(kind=dp), DIMENSION(2), INTENT(in) :: w
    COMPLEX(kind=dp), DIMENSION(2)             :: f
    
    !Calculate the components of the ODE
    f(1) = w(2)
    f(2) = 2*(V(x)-E)*w(1)
    
    RETURN
  END FUNCTION f

  
  FUNCTION V(x)
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


    !===============================DECLARATIONS================================!
    REAL(kind=dp) :: V
    REAL(kind=dp), INTENT(in) :: x
    
    !Triangle Potential
    IF (pot_form=='triangle') THEN
       IF (x .lt. a) THEN
          V=0.0_dp
       ELSE IF (x .ge. a .and. x .le. b) THEN
          V=((bias-height)*x)/(b-a) + height - ((bias-height)*a)/(b-a)
       ELSE IF (x .gt. b) THEN
          V=bias
       END IF

    !Parabolic triangle potential
    ELSE IF (pot_form=='parabolic_triangle') THEN
       IF (x .lt. a) THEN
          V=0.0_dp
       ELSE IF (x .ge. a .and. x .le. b) THEN
          V=para_cons*x**2 - 2.0_dp*b*para_cons*x - bias + para_cons*b**2
       ELSE IF (x .gt. b) THEN
          V=-bias
       END IF
    
    !Rectangle potential   
    ELSE IF (pot_form=='rectangle') THEN
       IF (x .lt. a) THEN
          V = 0.0_dp
       ELSE IF (x .ge. a .and. x .le. b) THEN
          V = (bias*x)/(b-a) + height - (bias*a)/(b-a)
       ELSE IF (x > b) THEN
          V = bias
       END IF
       
    !Double rectangle potential   
    ELSE IF (pot_form=='double_rectangle') THEN
       IF (x .lt. a) THEN
          V=0.0_dp
       ELSE IF (x .ge. a .and. x .le. b) THEN
          V=height
       ELSE IF (x .gt. b .and. x .lt. c) THEN
          V=0.0_dp
       ELSE IF (x .ge. c .and. x .le. d) THEN
          V=height
       ELSE IF (x .gt. d) THEN
          V=0.0_dp
       END IF
          
    !Trapezoid potential      
    ELSE IF (pot_form=='trapezoid') THEN
       IF (x .lt. a) THEN
          V=0.0_dp
       ELSE IF (x .ge. a .and. x .le. b) THEN
          V= height*(x-a)/(b-a)
       ELSE IF (x .gt. b .and. x .lt. c) THEN
          V=height
       ELSE IF (x .ge. c .and. x .le. d) THEN
          V= height + height*(x-c)/(c-d)
       ELSE IF (x .gt. d) THEN
          V=0.0_dp
       END IF

    !If the potential shape isn't recognised   
    ELSE
       STOP "Error finding potential for the input shape"
    END IF


    RETURN
  END FUNCTION V

  SUBROUTINE find_coeff
    IMPLICIT NONE
    !===========================================================================!
    ! This subroutine goes through the newly created output.dat file and gets   !
    ! the values of psi at x=0 and x=pi/2*k. Then, the values are used to       !
    ! calculate the reflection and transmission coefficients.                   !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! output.dat : data file created by runge_kutta subroutine                  !
    !---------------------------------------------------------------------------!                                                                          !
    ! Input:                                                                    !
    ! E          : Energy (Hartrees)                                            !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! T, R       : Transmission / Reflection coefficients                       !  
    !===========================================================================!


    !==============================DECLARATIONS=================================!            
    COMPLEX(kind=dp) :: psi_1,psi_2,At,Ar
    REAL(kind=dp)    :: x,Re,Im,k
    
    !Open the file to read the data from
    OPEN(15,file='output.dat',status='old',iostat=ierr)
    IF (ierr/=0) STOP "Error in opening output.dat"
    
    !Find psi at x = 0
    READ(15,*) x,Re,Im
    loop1 : DO
       IF (ABS(x) .gt. 1.0e-6_dp) THEN
          READ(15,*) x,Re,Im
       ELSE
          exit loop1
       END IF
    END DO loop1

    !Store the value
    psi_1 = CMPLX(Re,Im,kind=dp)
    
    !Close the file
    CLOSE(15,iostat=ierr)
    IF (ierr/=0) STOP "Error in close output.dat"

    !Open the file again to read the data from
    OPEN(15,file='output.dat',status='old',iostat=ierr)
    IF (ierr/=0) STOP "Error in opening output.dat"
    
    !Find psi at x = pi/2*k
    k = SQRT(2*E)
    READ(15,*) x,Re,Im
    loop2 : do
       IF (ABS(x) .gt. pi*0.5_dp/k) THEN
          READ(15,*) x,Re,Im
       ELSE
          exit loop2
       END IF
    END DO loop2
    
    !Store the value
    psi_2 = CMPLX(Re,Im,kind=dp)
    
    !Close the file
    CLOSE(15,iostat=ierr)
    IF (ierr/=0) STOP "Error in close output.dat"

    !Calculate R
    At   = (psi_1 - CMPLX(0.0_dp,1.0_dp,kind=dp)*psi_2)*0.5_dp
    Ar   = (CMPLX(0.0_dp,1.0_dp,kind=dp)*psi_2 + psi_1)*0.5_dp    
    refl = (ABS(Ar)**2)/(ABS(At)**2)

    RETURN
  END SUBROUTINE find_coeff
    
END MODULE ode_schrod
