##### ========================================
##### ===== CFSF FUNCTIONS AND CONSTANTS =====
##### ========================================

# the following submodule defines the (new) functions and constants
# that are used within the 'lhs' and 'function' arguments in create
# statements

# NOTE : (new) means that no direct alternative exists in Maple
# but this functionality can also be used to create aliasses
# of function names for use with query (e.g. uincgamma)

##### ----- SUBMODULE FUNCTIONS -----

functions := module()

  description "Functions and constants for use with the CFSF package";

  export
    # FUNCTIONS
    argument, pow, ln1p,
    Binet, Besseli, Besselj, Ein, WhittakerPsi,
    incBeta, regincBeta, lincgamma, uincgamma,
    cerf, erfcI, cnormalI, qhyper,
    normaldcdf, cnormaldcdf, gammacdf, cgammacdf, Mills,

    # CONSTANTS
    goldenratio, Gompertz, rabbit,

    # helper functions
    hankelsymb;

  options
    package;

##### ----- DEFINITIONS FOR FUNCTIONS -----

  argument := z -> `if`( type( z, zerotype ), undefined, `if`( type( z, symboltype ), 'functions:-argument(z)', :-argument( z ) ) );

  pow := (x,y) -> x^y;

  ln1p := z -> ln(1+z):

  Binet := z -> lnGAMMA(z) - (z-1/2) * ln(z) + z - ln(sqrt(2*Pi)):

  Besseli := (n,z) -> sqrt(Pi/(2*z)) * BesselI(n+1/2,z):

  Besselj := (n,z) -> sqrt(Pi/(2*z)) * BesselJ(n+1/2,z):

  Ein := z -> Ei(1,z) + :-gamma + ln(z):

  WhittakerPsi := (alpha,beta,z) -> z^((alpha+beta)/2-1) * exp(z/2) * WhittakerW(-(alpha+beta)/2,(beta-alpha)/2,z);

  incBeta := (z,a,b) -> ''int'( t^(a-1) * (1-t)^(b-1), t=0..z )';

  regincBeta := (z,a,b) -> incBeta(z,a,b) / Beta(a,b);

  lincgamma := (alpha,z) -> GAMMA(alpha) - GAMMA(alpha,z):

  uincgamma := (alpha,z) -> GAMMA(alpha, z):

  cerf := z -> exp(-z^2) * erfc(-I*z):

  erfcI := (alpha,z) -> erfc(alpha, z);

  qhyper := (a,b,c,z) -> 'qhyper(a,b,c,d)';  # TO BE ADDED (currently unevaluated)

  cnormalI := (k,x) -> 2^(k/2-1) * erfcI(k,x/sqrt(2));

  normaldcdf := x -> stats[statevalf,cdf,normald[0,1]](x);

  cnormaldcdf := x -> 1 - normaldcdf(x);

  gammacdf := (x,alpha,theta) -> ''stats[statevalf,cdf,gamma[alpha,theta]]'(x)';

  cgammacdf := (x,alpha,theta) -> 1 - gammacdf(x,alpha,theta);

  Mills := x -> sqrt(2*Pi) * exp(x^2/2) * cnormaldcdf(x);

##### ----- DEFINITIONS FOR CONSTANTS -----

  # the (new) constants are defined as functions without parameters

  goldenratio := () -> ( 1 + sqrt(5) ) / 2;

  Gompertz := () -> exp(1) * Ei(1,1);

  # NOTE : type( rabbit, constant ) does not evaluate to true;
  # to make it of type constant, we could add it to the global
  # 'constants' variable as well as define evalf/constant/rabbit
  # (likewise for goldenratio and Gompertz)

  rabbit := () -> sum( 2^( - floor( k * ( sqrt(5) + 1 ) / 2 ) ), k=1..infinity );

##### ----- definitions for helper functions -----

  # hankelsymb := (alpha,n) -> GAMMA(1/2+alpha+n) / ( n! * GAMMA(1/2+alpha-n) );
  hankelsymb := (alpha,n) -> (1/n!) * ''product'( ( alpha + t - 1/2 ) * ( alpha - t + 1/2 ), t=1..n )';

end module:
