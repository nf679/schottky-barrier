-------------------------WHAT THE PROGRAM DOES-----------------------------------------------------------------------

The program reads in a value of E but manually has input values of a starting energy of 0.015 with a number of
1000 iterations and a step size of 1e-5. It loops over this to produce data of E vs T.
There is a scaling parameter of E_00=1 added.


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
    energy        ::   Energy of the incident electron


---------------------------OUTPUT FILES-------------------------------------------------------------------------------

-output.dat is created after every run. It is overwritten if it already exists (as data is overwritten).
-Energy_v_coeff.dat is created after every run. It should be deleted before re-running the program (as data is appended).
