-------------------------WHAT THE PROGRAM DOES-----------------------------------------------------------------------

The program provides a data file of x, psi and psi' for every run.
It can optional provide a data file of V(x), x to show the effect of bias.
Otherwise, it takes in a bias range to loop over, calculates T then plots both
a file of bias vs T and Electric_field vs Current_Density. (I-V characteristics)


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
    energy        ::   Energy of the incident electrons
    bias_i        ::   Initial bias
    bias_f        ::   Final bias
    bias_n        ::   Number of iterations for bias calculation
    num_c         ::   Number of incident charge carriers (electrons)


---------------------------OUTPUT FILES-------------------------------------------------------------------------------

-output.dat is created after every run. It is overwritten if it already exists
-bias.dat is created after every run. It is overwritten if it already exists.
-current.dat is created after every run. It is overwritten if it already exists.


OPTIONAL:
-potential.dat can be made if uncommented. It is overwritten if it already exists.
