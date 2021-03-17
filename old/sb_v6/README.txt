-------------------------WHAT THE PROGRAM DOES-----------------------------------------------------------------------

The program asks the user for a start and end value of energy, as well as the number of iterations, and loops over 
this range to produce data for E vs transmission coefficient. It also provides a data file of Padovanni and 
Strattons model for E vs T, and a data file of x, psi and psi'
There is a parameter called sc_parameter which simply scales the P+S transmission coefficient.


----------------------------TO COMPILE:------------------------------------------------------------------------------

For normal run (using TISE):
  -Compile runge_kutta.f90 with ode_schrod.f90  


----------------------------INPUT FILES------------------------------------------------------------------------------

To run the program we require:

  -input1.dat
    start         ::   Position in space to start the program
    finish        ::   Position in space to finish the program
    n             ::   Number of iterations 

  -input2.dat
    shape         ::   The shape of the potential barrier
    height        ::   Height (energy/potential) of the potential barrier
    a b [c] [d]   ::   Start and end positions in space of the potential barrier [these are for the double rectangle]


---------------------------OUTPUT FILES-------------------------------------------------------------------------------

-output.dat is created after every run. It is overwritten if it already exists (as data is overwritten).
-modeltest.dat is created after every run. It should be deleted before re-running the program (as data is appended).
-Energy_v_coeff.dat is created after every run. It should be deleted before re-running the program (as data is appended).
