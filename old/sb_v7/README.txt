-------------------------WHAT THE PROGRAM DOES-----------------------------------------------------------------------

The program provides data files of x,psi,psi' and x,V(x) to illustrate the effect of an applied bias on the 
potential.

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
    bias          ::   Value of the applied bias (difference in potential on the rhs of the well)
    

---------------------------OUTPUT FILES-------------------------------------------------------------------------------

-output.dat is created after every run. It is overwritten if it already exists (as data is overwritten).
-potential.dat is created after every run. It is overwritten if it already exists (as data is overwritten).

