MODULE rk
  
  USE ode_schrod
  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE runge_kutta
    IMPLICIT NONE 
    !===========================================================================!
    ! This is an algorithm for the fourth order Runge Kutta method.             !
    ! This method is written so that it starts on b and moves backwards.        !
    !---------------------------------------------------------------------------!
    ! Requires:                                                                 !
    ! ode_schrod must be called before running this subroutine                  !
    !---------------------------------------------------------------------------!
    ! Input:                                                                    !
    ! input1.dat     : data file containing start, finish, n                    !
    !---------------------------------------------------------------------------!
    ! Output:                                                                   !
    ! output.dat     : data file containing x, psi, psi'                        !
    ! potential.dat  : data file containing x, V(x)                             !
    !===========================================================================!
    
    !===============================DECLARATIONS================================!
    INTEGER                        :: i,n,ierr
    REAL(kind=dp)                  :: x,h,start,finish,wave_n, x_temp
    COMPLEX(kind=dp), DIMENSION(2) :: w,k1,k2,k3,k4
    
    !Open the input file
    OPEN(10,file='input1.dat',status='old',iostat=ierr)
    IF (ierr/=0) STOP "Error in opening input1.dat"

    !Read the input parameters 
    READ(10,*) start,finish
    READ(10,*) n

    !Close the inpput file
    CLOSE(unit=10,iostat=ierr)
    IF (ierr/=0) STOP "Error in closing input1.dat"

    !Calculate starting psi and psi'
    wave_n = SQRT(2.0_dp*E)
    w(1)   = EXP(CMPLX(0.0_dp,wave_n,kind=dp)*finish)
    w(2)   = CMPLX(0.0_dp,wave_n,kind=dp)*w(1)        
    h      = (finish-start)/REAL(n,dp)
    x      = finish

    !Open a file to write the data to
    OPEN(15,file='output.dat',status='replace',iostat=ierr)
    IF (ierr/=0) STOP "Error in opening output.dat"

    !Write the data to the file
    WRITE(15,*) x,REAL(w(1)),AIMAG(w(1)),REAL(w(2)),AIMAG(w(2))
    
    !Runge kutta method
    DO i=1,n
       k1 = h*f(x,w)
       k2 = h*f(x-(h/2.0_dp),w-(k1/2.0_dp))
       k3 = h*f(x-(h/2.0_dp),w-(k2/2.0_dp))       
       k4 = h*f(x-h,w-k3)
       w  = w - (k1+2.0_dp*k2+2.0_dp*k3+k4)/6.0_dp
       x  = x - h
       !Write to the file
       WRITE(15,*) x,REAL(w(1)),AIMAG(w(1)),REAL(w(2)),AIMAG(w(2))
    END DO

    !Close the file the data was written to
    CLOSE(15,iostat=ierr)
    IF (ierr/=0) STOP "Error in closing output.dat"

    !Open a file to write the potential data to
    OPEN(35,file='potential.dat',status='replace',iostat=ierr)
    IF (ierr/=0) STOP "Error in opening potential.dat"

    !Set the starting x position to the start 
    x_temp = start

    !Loop over the x values to write v and x to the file
    DO i=1, n
      WRITE(35,*) x_temp, v(x_temp)
      x_temp = x_temp + h
    END DO
    
    !Close the potential data file 
    CLOSE(unit=35, iostat=ierr)
    IF (ierr/=0) STOP "Error in closing potential.dat"

    RETURN
  END SUBROUTINE runge_kutta
        
END MODULE rk


PROGRAM ode_solver
  USE rk  
  IMPLICIT NONE

  !=============================================================================!
  ! The main program which calls all of the subroutines                         !
  !=============================================================================!

  INTEGER :: i
  
  !Read the parameters for the bias and potential from the file
  CALL read_pot

  !Open a file to write the data to
  OPEN(55,file='bias.dat',status='replace',iostat=ierr)
  IF (ierr/=0) STOP "Error in opening bias.dat"
  
  E = E + bias 

  !Loop over the bias values
  DO i=1,bias_n

    !generate psi x data (depends on f(x) which depends on V(x) )
    CALL runge_kutta

    !calculate T, using psi x data
    CALL find_coeff

    !Write the stuff to the file
    WRITE(55,*) bias, 1.0_dp-refl

    !Increment the bias and energy
    bias = bias + bias_h

    E = E + bias_h

  END DO

  !Close the file the data was written to
  CLOSE(55,iostat=ierr)
  IF (ierr/=0) STOP "Error in closing bias.dat"
  
END PROGRAM ode_solver
