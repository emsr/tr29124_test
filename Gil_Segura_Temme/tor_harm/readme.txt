
       This package contains:
       
       *  DTORH1.F    This routine calculates Toroidal Harmonics for a given
                      order m, from the lowest degree (n=0) to a maximum
                      degree (n=N).
                       
       *  DTORH2.F    This routine calculates Toroidal Harmonics for orders
                      m=0,...,M and degrees n=0,...,N(m).

                       
       *  DTORH3.F    This routine calculates Toroidal Harmonics for orders
                      m=0,...,M and degrees n=0,...,NMAX, where NMAX is
                      the minimum value of the maximum degrees for each
                      order m.


       *  TESTPRO.F   Program test of the previous routines.
       
                      
       *  EXAMPLE.F   Program for the physical example mentioned in the
                      paper: Evaluation of the electrostatic field
                      due to a charged toroidal conductor at a potential
                      of the type: V=cos(m*phi).


       *  ROUT.F      Auxiliary routines referenced in DTORH1, DTORH2
                      and DTORH3. Then, they are needed also in TESTPRO.F
                      and EXAMPLE.F.

       *  TEST.OUT    Output file of the test program.           
       
       
