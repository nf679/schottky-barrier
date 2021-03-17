To compile:

   For TEST run:
       -Compile ode_test.f90 (module)

       -In runge_kutta.f90 chage third line to "use ode_test"
       Comment out lines 79 & 81  ("call read_pot" and "call find_coeff")
       
       -Compile runge_kutta.f90 with ode_test.o

   For normal run (using TISE):
       -Compile ode_schrod.f90 (module)

       -In runge_kutta.f90 chage third line to "use ode_schrod"
       Lines 79 & 81  ("call read_pot" and "call find_coeff") should NOT
       be commented out

       -Compile runge_kutta.f90 with ode_schrod.o

Input Files:

   -All runs require input1.dat

   -normal runs (using TISE from ode_schrod.f90) require input2.dat

Output:

   -Reflection Coefficient is printed to screen

   -Output.dat is created after every run. Note, it is overwritten if it already
    exists
