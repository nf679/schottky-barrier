-------------------------WHAT THE PROGRAM DOES-----------------------------------------------------------------------

The program has a number of output options. For every option the program creates a file with x, psi data in.

----------------------------TO COMPILE:------------------------------------------------------------------------------

For normal run (using TISE):
  -Compile runge_kutta.f90 with ode_schrod.f90  

---------------------------OUTPUT OPTIONS-------------------------------------------------------------------------------

1: position (x) vs wavefuction (psi)
2: position (x) vs first derivative of wavefunction (psi')
3: Energy of incident electron (E) vs Transmission Coefficient (T)
4: Energy of incident electron (E) vs Reflection Coefficient (R) 
5: Position (x) vs Potential (v)


----------------------------INPUT FILES------------------------------------------------------------------------------

To run the program we require:

  -input1.dat
    start         ::   Position in space to start the program
    finish        ::   Position in space to finish the program
    n             ::   Number of iterations 

  -input2.dat
    output_option ::   Depending on what we want data for
    shape         ::   The shape of the potential barrier
    height        ::   Height (energy/potential) of the potential barrier
    a b [c] [d]   ::   Start and end positions in space of the potential barrier [these are for the double rectangle]
    energy        ::   Energy of the incident electron (If option 1,2, or 5 is chosen)
    e_start       ::   Starting energy (if option 3 or 4 is chosen)
    e_finish      ::   Final energy (if option 3 or 4 is chosen)
    e_iterations  ::   Number of iterations for the energy calculation (if option 3 or 4 is chosen)


---------------------------OUTPUT FILES-------------------------------------------------------------------------------

-output1.dat is created after every run. It is overwritten if it already exists (as data is overwritten).

OPTIONAL

-output2.dat overwritten if it already exists
-output3.dat appended so delete if rerunning program  
-output4.dat appended so delete if rerunning program 
-output5.dat overwritten if it already exists