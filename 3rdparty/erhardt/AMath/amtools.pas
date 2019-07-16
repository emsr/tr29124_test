unit AMTools;

{Accurate Math Tools: Root finding, function minimization, numerical
                      integration, convergence acceleration,
                      quadratic/cubic equations}

interface

{$i STD.INC}

{$ifdef BIT16}
{$N+,F+}
{$endif}

uses
  amath;

(*************************************************************************

 DESCRIPTION   :  Accurate Math Tools: Root finding, function minimization, numerical
                  integration, convergence acceleration, quadratic/cubic equations

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    : References used in this unit, main index in amath_info.txt/references

                 [13] W.H. Press et al, Numerical Recipes in C, 2nd ed., Cambridge, 1992,
                      http://www.nrbook.com/a/bookcpdf.html
                 [15] W. Kahan, On the Cost of Floating-Point Computation
                      Without Extra-Precise Arithmetic, 2004.
                      http://www.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf
                 [16] G.E. Forsythe, How do you solve a quadratic equation?
                      Stanford University Technical Report no. CS40, 1966
                 [17] P.H. Sterbenz, Floating-Point Computation, 1974, Chap.9.3
                      Carefully written programs/quadratic equation, p.246ff
                 [18] W.Y. Sit, Quadratic Programming? 1997.
                      http://www.mmrc.iss.ac.cn/ascm/ascm03/sample.pdf
                 [28] R.P. Brent, Algorithms for Minimization without Derivatives,
                      Englewood Cliffs, 1973. Scanned copy available from the author's
                      site: http://maths-people.anu.edu.au/~brent/pub/pub011.html
                 [29] G.E. Forsythe, M.A. Malcolm, C.B. Moler, Computer Methods for
                      Mathematical Computations, Englewood Cliffs, 1977.
                      Fortran code from http://www.netlib.org/fmm/
                 [34] R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D. Kahaner,
                      QUADPACK: A subroutine package for automatic integration (1983).
                      Public domain Fortran source: http://www.netlib.org/quadpack/
                 [38] T. Ooura's Fortran and C source code for automatic quadrature
                      using Double Exponential transformation; available from
                      http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html
                 [56] E.J. Weniger, Nonlinear sequence transformations for the acceleration
                      of convergence and the summation of divergent series, 1989,
                      Comput. Phys. Rep. 10, pp.189-371; available as
                      https://arxiv.org/pdf/math/0306302v1.pdf
                 [57] T. Fessler, W.F. Ford, D.A. Smith, HURRY: An Acceleration Algorithm for
                      Scalar Sequences and Series, ACM TOMS, Vol.9, No.3, 1983, pp.346-354.
                      Fortran source available from http://netlib.org/toms/602
                 [60] W. Kahan, To Solve a Real Cubic Equation (Lecture Notes
                      for a Numerical Analysis Course), 1986. Available as
                      http://www.dtic.mil/dtic/tr/fulltext/u2/a206859.pdf

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     01.10.10  W.Ehrhardt  Initial BP7 version of zbrent
 0.11     02.10.10  we          Error codes for zbrent
 0.12     02.10.10  we          zeroin
 0.13     03.10.10  we          First version of localmin
 0.14     04.10.10  we          Simplified version mbrent and fmin
 0.15     10.03.11  we          quanc8
 0.16     11.03.11  we          quanc8: make tolerances non-negative
 0.17     17.03.11  we          qk21
 0.18     17.03.11  we          qagse, qpsrt, qelg
 0.19     18.03.11  we          qk15i, qagie
 0.20     18.03.11  we          am_getmem, am_freemem, qags, qagi
 0.21     19.03.11  we          quagk
 0.22     28.03.11  we          Fix Delphi 2 optimization bug
 0.23     31.03.11  we          rewrite qpsrt
 0.24     28.09.11  we          qawce/qawc
 0.25     05.10.11  we          intdei/intdeo
 0.26     15.10.11  we          intde
 0.27     11.08.13  we          Change type of lsx to word (for 16-bit limit > 3200)
 0.28     22.08.13  we          Levin u-transformation from whiz1
 0.29     23.08.13  we          Wynn epsilon
 0.30     24.08.13  we          Renamed whiz1 to levinu1 (no sequence)
 0.31     24.08.13  we          wynneps1 with sum parameter
 0.32     25.08.13  we          levinu1/wynneps1 for $ifndef CONST, param descriptions
 0.33     23.09.13  we          cubsolve
 0.34     15.04.14  we          am_getmem/freemem replaced with memh.malloc/mfree
 0.35     05.10.16  we          zbrenty for zeros of f(x)-y
 0.36     29.06.17  we          Code removed for TP5-TP6, TPW1-D1
 0.37     30.06.17  we          PolyRoots functions
 0.38     09.07.17  we          PolyRoots polishes complex roots
 0.39     28.07.17  we          Removed PolyRootSort, PolyRootPolish from interface
 0.40     02.12.17  we          Suppress warnings: Local variable does not seem to be initialized
 0.41     08.01.18  we          intde_p and intdei_p
 0.42     02.02.18  we          Removed ier=2 for intde* routines
 0.43     03.02.18  we          PolyRootsOA (open array version)
 0.44     12.02.18  we          Avoid non-sense warning for FPC304 -O4
 0.45     19.04.18  we          Fixes for FPC311
 0.46     26.10.18  we          Ridders' root finding method zridders

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2009-2018 Wolfgang Ehrhardt

 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from
 the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software in
    a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution.
----------------------------------------------------------------------------*)

{#Z+}
{---------------------------------------------------------------------------}
{----------------------- Function minimization -----------------------------}
{---------------------------------------------------------------------------}
{#Z-}
procedure localmin(f: TFuncX; a,b,eps,t: extended; var x,fx: extended; var ic: integer);
  {-Brent's algorithm (with guaranteed convergence) for finding a local    }
  { minimum of the function f in the interval (a,b). x is the approximate  }
  { minimum abscissa, fx=f(x). eps and t define a tolerance tol =eps*|x|+t.}
  { f is never evaluated for 2 points closer together than tol. eps shall  }
  { not be < 2*eps_x, preferably not smaller than sqrt(eps_x). ic is the   }
  { iteration count, -1 if a=b, 0 if max. count = 5000 exceeded.}
  { The algorithm combines golden section search and successive parabolic  }
  { interpolation using only function (not derivative) evaluations.        }

procedure mbrent(f: TFuncX; a,b,t: extended; var x,fx: extended; var ic: integer);
  {-Find a local minimum of the function f in the interval (a,b). }
  { x is the approximate minimum abscissa, fx = f(x). Simplified  }
  { version of procedure localmin with fixed eps=0.5*sqrt(eps_x). }
  { ic is the iteration count, -1 if a=b, 0 if max.=5000 exceeded.}

function fmin(f: TFuncX; a,b,tol: extended): extended;
  {-Find a local minimum of the function f in the interval (a,b) and return }
  { the approximate minimum abscissa. fmin is a simple shell version for the}
  { procedure localmin using eps = 0.5*sqrt(eps_x) and t = tol/3.}

{#Z+}
{---------------------------------------------------------------------------}
{--------------------------- Root finding ----------------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function zbrenty(f: TFuncX; y,a,b,t: extended; var ic,err: integer): extended;
  {-Brent/Dekker algorithm with guaranteed convergence for finding a zero   }
  { x of the function f(x)-y in the interval [a,b] to within a tolerance of }
  { 6*eps_x*|x|+2*t, where t is a positive tolerance, i.e. it solves f(x)=y.}
  { Assumes that f(a)-y and f(b)-y have different signs. ic is the iteration}
  { count; err is an error code (0: no error, -1: if f(a)-y and f(b)-y have }
  { the same sign, -2: max. iteration count exceeded). The algorithm is     }
  { based on a combination of successive interpolations and bisection.      }

function zbrent(f: TFuncX; a,b,t: extended; var ic,err: integer): extended;
  {-Brent/Dekker algorithm with guaranteed convergence for finding a zero  }
  { of a function: Return a zero x of the function f in the interval [a,b],}
  { function zbrenty with y=0.0 is used.                                   }

function zeroin(f: TFuncX; a,b,t: extended): extended;
  {-Return a zero x of the function f in the interval [a, b] to within a  }
  { tolerance 6*eps_x*|x| + 2*t, where t is a positive tolerance, assumes }
  { that f(a) and f(b) have different signs. Simplified version of zbrent.}

function zridders(f: TFuncX; a,b,t: extended): extended;
  {-Find a zero x of f with Ridders' method in the interval [a,b] with}
  { a tolerance 2*eps_x*|x| + 0.5*t, returns NaN if f(a)*f(b) > 0.    }

{#Z+}
{---------------------------------------------------------------------------}
{--------------------- Numerical integration -------------------------------}
{---------------------------------------------------------------------------}
{#Z-}
procedure quanc8(fun: TFuncX; a,b,abserr,relerr: extended;
                 var result, errest, flag: extended; var nofun: longint);
  {-Estimate the integral of fun(x) from a to b to a user provided tolerance.}
  { Pascal translation of the Fortran subroutine by Forsythe, Malcolm, Moler.}

  {#F}
  {An automatic adaptive routine based on the 8-panel Newton-Cotes rule: }
  {Input:                                                                }
  {  fun     The name of the integrand function subprogram fun(x).       }
  {  a       The lower limit of integration.                             }
  {  b       The upper limit of integration, b may be less than a.       }
  {  relerr  A  relative error tolerance,  should be non-negative.       }
  {  abserr  An absolute error tolerance,  should be non-negative.       }
  {Output:                                                               }
  {  result  An approximation to the integral hopefully satisfying the   }
  {          least stringent of the two error tolerances.                }
  {  errest  An estimate of the magnitude of the actual error.           }
  {  nofun   The number of function values used in calculation of result.}
  {  flag    A reliability indicator: if flag is zero, then result       }
  {          probably satisfies the error tolerance. If flag is          }
  {          xxx.yyy , then xxx = the number of intervals which have.    }
  {          not converged and 0.yyy = the fraction of the interval      }
  {          left to do when the limit on nofun was approached.          }
  {#F}


const
  QMAXLIM = 5000;    {Maximum limit of subintervals for adaptive quadrature}

var
  DefLimit: integer; {Default limit = 500}

type
  TQXLimArr = array[1..QMAXLIM] of extended;  {Extended vector for subinterval lists}
  TQILimArr = array[1..QMAXLIM] of integer;   {Integer  vector for subinterval lists}

procedure quagk(f: TFuncX; a,b,epsabs: extended; var result,abserr: extended; var ier: integer);
  {-Global adaptive quadrature of f over (a,b) based on Gauss-Kronrod rules}
  { for the subintervals, with acceleration by Wynn's epsilon algorithm.}

  {#F}
  { Simplified user interface to procedure qags and qagi.}
  { f      - function defining the integrand                       }
  { a      - lower limit of integration, may be infinite           }
  { b      - upper limit of integration, may be infinite           }
  { epsabs - absolute accuracy requested                           }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { ier    - return code (success=0, error>0 see below)            }

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine the estimates for integral  }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges. If possible, }
  {        an appropriate special-purpose integrator should be used, which }
  {        is designed for handling the type of difficulty involved.       }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved. The error may be       }
  {        under-estimated.                                                }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 4 The algorithm does not converge. Roundoff error is detected in  }
  {        the extrapolation table. It is presumed that the requested      }
  {        tolerance cannot be achieved, and that the returned result is   }
  {        the best which can be obtained.                                 }
  {                                                                        }
  {    = 5 The integral is probably divergent, or slowly convergent. It    }
  {        must be noted that divergence can occur with any other value of }
  {        ier.                                                            }
  {                                                                        }
  {    = 6 The input is invalid, because epsabs <= 0 and epsrel < 50*eps_x.}
  {        result, abserr, last are set to zero.                           }
  {                                                                        }
  {    = 7 The input is invalid, limit < 0 or limit > QMAXLIM.             }
  {        result, abserr, last are set to zero.                           }
  {                                                                        }
  {    = 8 Dynamic list vectors cannot be allocated.                       }
  {        result, abserr, last are set to zero.                           }
  {                                                                        }
  {    = 9 At least one limit a or b is NaN, or infinite a=b.              }
  {        result, abserr, last are set to NaN_x.                          }
  {#F}

procedure qags(f: TFuncX; a, b, epsabs, epsrel: extended; limit: integer;
               var result, abserr: extended; var neval: longint; var ier: integer);
  {-Global adaptive quadrature of f over (a,b) based on 21-point Gauss-Kronrod}
  { rule for the subintervals, with acceleration by Wynn's epsilon algorithm.}
  { Simplified user interface to procedure qagse}

  {#F}
  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the no. of subintervals, 0: use DefLimit}
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  {#F}

procedure qagi(f: TFuncX; bound: extended;inf: integer; epsabs, epsrel: extended; limit: integer;
               var result, abserr: extended; var neval: longint; var ier: integer);
  {-Global adaptive quadrature of f over an infinite interval based on trans- }
  { formed 15-point Gauss-Kronrod rule for the subintervals, with acceleration}
  { by Wynn's epsilon algorithm. Simplified user interface to procedure qagie.}

  {#F}
  { f      - function defining the integrand                       }
  { bound  - finite bound of integration ran                       }
  { inf    - indicating the kind of integration range involved     }
  {              inf =  1 corresponds to  (bound, +infinity)       }
  {              inf = -1             to  (-infinity, bound)       }
  {              inf =  2             to  (-infinity, +infinity)   }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the no. of subintervals, 0: use DefLimit}
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  {#F}

procedure qagse(f: TFuncX; a, b, epsabs, epsrel: extended;
                limit: integer; var result, abserr: extended;
                var neval: longint; var ier: integer;
                var alist, blist, elist, rlist: TQXLimArr;
                var iord: TQILimArr; var last: integer);
  {-Global adaptive quadrature of f over (a,b) based on 21-point Gauss-Kronrod}
  { rule for the subintervals, with acceleration by Wynn's epsilon algorithm.}

  {#F}
  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the number of subintervals              }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  { alist  - left end points of the subintervals                   }
  { blist  - right end points of the subintervals                  }
  { rlist  - integral approximations on the subintervals           }
  { elist  - absolute error estimates on the subintervals          }
  { iord   - first elements contain pointers to the error estimates}
  { last   - number of subintervals actually produced              }
  {#F}

procedure qagie(f: TFuncX; bound: extended; inf: integer;epsabs, epsrel: extended; limit: integer;
                var result, abserr: extended;
                var neval: longint;
                var ier: integer;
                var alist, blist, elist, rlist: TQXLimArr;
                var iord: TQILimArr; var last: integer);
  {-Global adaptive quadrature of f over an infinite interval based on}
  { transformed 15-point Gauss-Kronrod rule for the subintervals, with}
  { acceleration by Wynn's epsilon algorithm.}

  {#F}
  { f      - function defining the integrand                       }
  { bound  - finite bound of integration ran                       }
  { inf    - indicating the kind of integration range involved     }
  {              inf =  1 corresponds to  (bound, +infinity)       }
  {              inf = -1             to  (-infinity, bound)       }
  {              inf =  2             to  (-infinity, +infinity)   }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the number of subintervals, must be > 1 }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  { alist  - left end points of the subintervals                   }
  { blist  - right end points of the subintervals                  }
  { rlist  - integral approximations on the subintervals           }
  { elist  - absolute error estimates on the subintervals          }
  { iord   - first elements contain pointers to the error estimates}
  { last   - number of subintervals actually produced              }
  {#F}

procedure qk21(f: TFuncX; a,b: extended; var result,abserr,resabs,resasc: extended);
  {-Integrate function f over (a,b) with 21-point Gauss-Kronrod rule}

  {#F}
  { a      - lower limit of integration                     }
  { b      - upper limit of integration                     }
  { result - approximation to the integral I of f over (a,b)}
  { abserr - estimate of the modulus of the absolute error  }
  { resabs - approximation to the integral |f| over (a,b)   }
  { resasc - approximation integral |f-I/(b-a)| over (a,b)  }
  {#F}

procedure qawce(f: TFuncX; a, b, c, epsabs, epsrel: extended;
                limit: integer; var result, abserr: extended;
                var neval: longint; var ier: integer;
                var alist, blist, elist, rlist: TQXLimArr;
                var iord: TQILimArr; var last: integer);
  {-Adaptive quadrature of the function f(x)/(x-c) over the finite interval}
  { (a,b) with the singularity at c with c<>a, c<>b. The routine calculates}
  { an approximation result to the Cauchy principal value. Parameters:}

  {#F}
  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { c      - singularity, ier=6 if c=a or c=b                      }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the number of subintervals              }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  { alist  - left end points of the subintervals                   }
  { blist  - right end points of the subintervals                  }
  { rlist  - integral approximations on the subintervals           }
  { elist  - absolute error estimates on the subintervals          }
  { iord   - first elements contain pointers to the error estimates}
  { last   - number of subintervals actually produced              }
  {#F}

procedure qawc(f: TFuncX; a, b, c, epsabs, epsrel: extended;
               limit: integer; var result, abserr: extended;
               var neval: longint; var ier: integer);
  {-Adaptive quadrature of the function f(x)/(x-c) over the finite interval}
  { (a,b) with the singularity at c with c<>a, c<>b. The routine calculates}
  { an approximation result to the Cauchy principal value.  Simplified user}
  { interface to procedure qawce. Parameters:}

  {#F}
  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { c      - singularity, ier=6 if c=a or c=b                      }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the no. of subintervals, 0: use DefLimit}
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  {#F}


procedure intde(f: TFuncX; a, b, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }

  {#F}
  { f      - integrand f(x), must be analytic over (a,b)   }
  { a      - lower limit of integration                    }
  { b      - upper limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}

{---------------------------------------------------------------------------}
procedure intde_p(f: TFuncXP; p: pointer; a, b, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }

  {#F}
  { f      - integrand f(x), must be analytic over (a,b)   }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { b      - upper limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}

procedure intdei(f: TFuncX; a, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { a      - lower limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}

procedure intdei_p(f: TFuncXP; p: pointer; a, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}

procedure intdeo(f: TFuncX; a, omega, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has an oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { a      - lower limit of integration                    }
  { omega  - frequency of oscillation, i.e. f has the form }
  {          f(x) = g(x)*sin(omega*x + theta) as x -> INF. }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32, eps>1, or omega=0   }
  {          3: max. iterations, result/abserr have values }
  {#F}


{#Z+}
{---------------------------------------------------------------------------}
{----------------------- Convergence acceleration --------------------------}
{---------------------------------------------------------------------------}
{#Z-}
{---------------------------------------------------------------------------}
procedure levinu1(an: extended; n,nmax: integer;
                  var qnum, qden: array of extended;
                  var extsum, sn: extended; var ierr: integer);
  {-Perform one step of the Levin u-transformation for sum(an,n=0..). }
  { The driver routine has to analyze the convergence of the process. }
  { DO NOT CHANGE the var parameters between successive calls!        }
  {#F}
  { an     - the n'th term of the sum                                 }
  { n      - the index of current term for this step                  }
  { nmax   - maximum valid value for n and index the a arrays         }
  { qnum   - backward diagonal of numerator array,   at least 0..nmax }
  { qden   - backward diagonal of denominator array, at least 0..nmax }
  { extsum - extrapolated value of the sum                            }
  { sn     - the n'th partial sum, will be set to an if n=0           }
  { ierr   - error if <> 0 (currently for n < 0,  n > nmax)           }
  {#F}

procedure wynneps1(an: extended; n,nmax: integer; sum: boolean;
                   var e: array of extended;
                   var extlim, sn: extended; var ierr: integer);
  {-Perform one step of Wynn's epsilon algorithm for the sequence  }
  { an, n=0.. or the sum(an,n=0..). Note: The driver routine has to}
  { analyze the convergence of the whole process. DO NOT CHANGE the}
  { var parameters between successive calls!                       }
  {#F}
  { an     - the n'th term ot the sequence or sum                  }
  { n      - the index of current term for this step               }
  { nmax   - maximum valid value for n and index the e array       }
  { sum    - if true then a sum is processed, otherwise a sequence }
  { e      - working array for epsilon table, at least 0..nmax     }
  { extlim - extrapolated value of the sequence or sum             }
  { sn     - the n'th partial sum, will be set to an if n=0        }
  { ierr   - error if <> 0 (currently for n < 0,  n > nmax)        }
  {#F}


{#Z+}
{---------------------------------------------------------------------------}
{----------------- Quadratic, cubic, polynomial equations ------------------}
{---------------------------------------------------------------------------}
{#Z-}

const
  MaxDeg = 20; {Maximum degree of polynomials for PolyRoots}

type
  TPolyVec = array[0..MaxDeg] of double; {Coefficients/roots of a polynom}

function squad(a,b,c: double; var x1,y1,x2,y2: double): integer;
  {-Solve the quadratic equation a*x^2 + b*x + c = 0. Result is the number}
  { of different solutions: 0 (if a=b=0), 1 (x1), or 2 (x1,x2). If the}
  { result is = -2, x1+i*y1 and x2+i*y2 are the two complex solutions.}
  { No precautions against over/underflow, NAN/INF coefficients return 0.}

function squadx(a,b,c: double; var x1,y1,x2,y2: double): integer;
  {-Solve the quadratic equation a*x^2 + b*x + c = 0. Result is the number}
  { of different solutions: 0 (if a=b=0 or INF/NAN), 1 (x1), or 2 (x1,x2).}
  { If the result is = -2, then x1 + i*y1 and x2 + i*y2 are the two complex}
  { solutions. Uses scaling by powers of two to minimize over/underflows.}

procedure cubsolve(a,b,c,d: double; var x,x1,y1,x2,y2: double);
  {-Solve the cubic equation ax^3 + bx^2 + cx + d = 0: compute a real root}
  { x (may be INF if a~0) and two complex zeros x1 + i*y1, x2 + i*y2 where}
  { y2 = -y1 may be zero, i.e. there are three reel roots.}

procedure PolyRoots(const p: TPolyVec; n: integer; var x,y: TPolyVec; var ierr: integer);
  {-Compute the n (complex) roots x[k] + i*y[k] of the polynomial  }
  { p(z) = p[0] + p[1]*z + ... p[n]*z^n, n>0, p[n]<>0. Real roots  }
  { have y[k]=0. ierr is an error code 0: OK, -1: n<1, -2: p[n]=0, }
  { k: iteration count exceeded while the kth root is being sought.}
  { n must be <= MaxDeg (currently = 20); the cases n=1 and 2 are  }
  { handled separately: n=1 directly and n=2 with squad,           }

  { PolyRoots uses a companion matrix method, balancing, and the   }
  { QR algorithm for the eigenvalues of an upper Hessenberg matrix.}

procedure PolyRootsA(const p: TPolyVec; n: integer; var x,y: TPolyVec; var ierr: integer);
  {-Perform PolyRoots, then improve and sort the roots}

procedure PolyRootsOA(const p: array of double; n: integer; var x,y: TPolyVec; var ierr: integer);
  {-Open array version, n <= high(p): perform PolyRoots, then improve and sort the roots}

{#Z+}
{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}
procedure qpsrt(limit: integer; var last, maxerr: integer; var ermax: extended;
                var elist: TQXLimArr; var iord: TQILimArr; var nrmax: integer);
  {-Maintain the descending ordering in the list of the local error estimates}

  { limit  - maximum number of error estimates the list can contain}
  { last   - number of error estimates currently in the list}
  { maxerr - maxerr points to the nrmax-th largest error estimate currently in the list}
  { ermax  - nrmax-th largest error estimate ermax = elist[maxerr]}
  { elist  - vector containing the error estimates}
  { iord   - the first k elements contain pointers to the error estimates}
  { nrmax  - maxerr = iord(nrmax)}

{#Z-}


implementation

uses
  memh;

{---------------------------------------------------------------------------}
{----------------------- Function minimization -----------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure localmin(f: TFuncX; a,b,eps,t: extended; var x,fx: extended; var ic: integer);
  {-Brent's algorithm (with guaranteed convergence) for finding a local    }
  { minimum of the function f in the interval (a,b). x is the approximate  }
  { minimum abscissa, fx=f(x). eps and t define a tolerance tol =eps*|x|+t.}
  { f is never evaluated for 2 points closer together than tol. eps shall  }
  { not be < 2*eps_x, preferably not smaller than sqrt(eps_x). ic is the   }
  { iteration count, -1 if a=b, 0 if max. count = 5000 exceeded.           }

  { The algorithm combines golden section search and successive parabolic  }
  { interpolation using only function (not derivative) evaluations.        }
const
  c = 0.381966011250105152; {(3 - sqrt(5))/2}
  ITMax = 5000;
var
  d,e,m,p,q,r,tol,t2,u,v,w: extended;
  fu,fv,fw: extended;
  it: integer;
begin

  {Ref: Brent [28] ch.5, [29] ch.8; see also [13] ch.10.2. The convergence}
  {in a reasonable number of function evaluations is guaranteed, for a C_2}
  {function with positive curvature at the minimum it will be superlinear }
  {(provided that the minimum is at an interior point of the interval).   }

  {Initialization}
  x  := a + c*(b - a);
  fx := f(x);
  if a=b then begin
    {Return with error indicator and x=a=b, fx = f(x)}
    ic := -1;
    exit;
  end;
  v  := x;
  w  := x;
  fv := fx;
  fw := fx;
  e  := 0.0;
  d  := 0.0;

  {Force t >= 0}
  t := abs(t);
  {eps shall not be < 2*eps_x, preferably not smaller than sqrt(eps_x)}
  if eps < 2.0*eps_x then eps := 2.0*eps_x;

  for it:=1 to ITMax do begin
    ic  := it;
    m   := 0.5*(a + b);
    tol := eps*abs(x) + t;
    t2  := 2.0*tol;

    {Check stopping criterion}
    if abs(x-m) <= t2 - 0.5*(b-a) then exit;

    if abs(e) <= tol then begin
      {A "golden-section" step}
      if x < m then e := b-x else e := a-x;
      d := c*e;
    end
    else begin
      {Fit parabola}
      r := (x-w)*(fx-fv);
      q := (x-v)*(fx-fw);
      p := (x-v)*q - (x-w)*r;
      q := 2.0*(q-r);
      if q>0.0 then p := -p else q := -q;
      r := e;
      e := d;
      {Note that there is a typo in Brent's Algol function (FORTRAN is OK): }
      {his condition p<q*(a-x) should be p>q*(a-x). This is the requirement,}
      {that the interpolated new u = x + p/q must be in the interval (a,b). }
      if (abs(p) < abs(0.5*q*r)) and (p > q*(a-x)) and (p < q*(b-x)) then begin
        {A "parabolic interpolation" step}
        d := p/q;
        u := x + d;
        {f must not be evaluated too close to a or b}
        if (u-a < t2) or (b-u < t2) then begin
          if x < m then d := tol else d := -tol;
        end;
      end
      else begin
        {A "golden-section" step}
        if x < m then e := b-x else e := a-x;
        d := c*e;
      end;
    end;

    {f must not be evaluated too close to x}
    if abs(d) > tol then u := x + d
    else if d > 0.0 then u := x + tol
    else u := x - tol;
    fu := f(u);

    {update a, b, v, w, and x}
    if fu <= fx then begin
      if u < x then b := x else a := x;
      v  := w;
      fv := fw;
      w  := x;
      fw := fx;
      x  := u;
      fx := fu;
    end
    else begin
      if u < x then a := u else b := u;
      if (fu <= fw) or (w=x) then begin
        v  := w;
        fv := fw;
        w  := u;
        fw := fu;
      end
      else if (fu <= fv) or (v=x) or (v=w) then begin
        v  := u;
        fv := fu;
      end;
    end;
  end;
  {Indicate no convergence}
  ic := 0;
end;


{---------------------------------------------------------------------------}
procedure mbrent(f: TFuncX; a,b,t: extended; var x,fx: extended; var ic: integer);
  {-Find a local minimum of the function f in the interval (a,b). }
  { x is the approximate minimum abscissa, fx = f(x). Simplified  }
  { version of procedure localmin with fixed eps=0.5*sqrt(eps_x). }
  { ic is the iteration count, -1 if a=b, 0 if max.=5000 exceeded.}
begin
  localmin(f,a,b,0.5*sqrt(eps_x),t,x,fx,ic);
end;


{---------------------------------------------------------------------------}
function fmin(f: TFuncX; a,b,tol: extended): extended;
  {-Find a local minimum of the function f in the interval (a,b) and return }
  { the approximate minimum abscissa. fmin is a simple shell version for the}
  { procedure localmin using eps = 0.5*sqrt(eps_x) and t = tol/3.}
var
  x,fx: extended;
  ic: integer;
begin
  localmin(f,a,b,0.5*sqrt(eps_x),tol/3.0,x,fx,ic);
  fmin := x;
end;


{---------------------------------------------------------------------------}
{--------------------------- Root finding ----------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function zbrenty(f: TFuncX; y,a,b,t: extended; var ic,err: integer): extended;
  {-Brent/Dekker algorithm with guaranteed convergence for finding a zero   }
  { x of the function f(x)-y in the interval [a,b] to within a tolerance of }
  { 6*eps_x*|x|+2*t, where t is a positive tolerance, i.e. it solves f(x)=y.}
  { Assumes that f(a)-y and f(b)-y have different signs. ic is the iteration}
  { count; err is an error code (0: no error, -1: if f(a)-y and f(b)-y have }
  { the same sign, -2: max. iteration count exceeded). The algorithm is     }
  { based on a combination of successive interpolations and bisection.      }

const
  ITMax = 5000;
var
  c,d,e,fa,fb,fc,tol,m,p,q,r,s: extended;
  it: integer;
begin

  {Ref: Brent [28] ch.3/4, [29] ch.7; see also [13] ch.9.3. This algorithm }
  {uses only function (not derivative) evaluations and never converges much}
  {more slowly than bisection. In practice the convergence is usually much }
  {faster than for bisection. If x is a simple zero of a C_1 function then }
  {superlinear convergence is finally achieved in many cases.              }

  {Force t >= 0}
  t := abs(t);

  {initialization}
  ic := 0;
  err:= 0;

  fa := f(a)-y;
  if fa=0.0 then begin
    zbrenty := a;
    exit;
  end;
  fb := f(b)-y;
  if fb=0.0 then begin
    zbrenty := b;
    exit;
  end;

  if (fa>0.0)=(fb>0.0) then begin
    zbrenty := b;
    err := -1;
    exit;
  end;

  {In the Brent/FMM code this is inside the loop and may be skipped via goto}
  {from the end. My Pascal code uses two copies of the lines to avoid goto.}
  c  := a;
  fc := fa;
  d  := b - a;
  e  := d;

  for it:=1 to ITMax do begin
    ic := it;
    if abs(fc) < abs(fb) then begin
      a  := b;
      b  := c;
      c  := a;
      fa := fb;
      fb := fc;
      fc := fa;
    end;
    {b and c bracket the zero, b is the latest iterate and the closest}
    {approximation, |f(b)| <= |f(c)|, and a is the previous iterate.}

    {convergence test, tol is >= 0}
    tol := 2.0*eps_x*abs(b) + 0.5*t;
    m   := 0.5*(c - b);
    if (abs(m) <= tol) or (fb = 0.0) then begin
      zbrenty := b;
      exit;
    end;

    {test if bisection is necessary}
    if (abs(e) < tol) or (abs(fa) <= abs(fb)) then begin
      {bisection}
      d := m;
      e := d;
    end
    else begin
      s := fb/fa;
      if a=c then begin
        {linear interpolation}
        p := 2.0*m*s;
        q := 1.0 - s;
      end
      else begin
        {inverse quadratic interpolation}
        q := fa/fc;
        r := fb/fc;
        p := s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
        q := (q - 1.0)*(r - 1.0)*(s - 1.0);
      end;

      {adjust signs}
      if p > 0.0 then q := -q else p := -p;

      {test if interpolation is acceptable}
      if (2.0*p >= 3.0*m*q - abs(tol*q)) or (p >= abs(0.5*e*q)) then begin
        {bisection}
        d := m;
        e := d;
      end
      else begin
        {use interpolation}
        e := d;
        d := p/q;
      end;
    end;

    {complete step, update previous (a,fa) and calculate new b and f(b)}
    a  := b;
    fa := fb;
    if abs(d) > tol then b := b + d
    else if m > 0.0 then b := b + tol
    else b := b - tol;
    fb := f(b)-y;
    if (fb>0.0)=(fc>0.0) then begin
      c  := a;
      fc := fa;
      d  := b - a;
      e  := d;
    end;
  end;
  {No convergence after ITMax steps}
  err := -2;
  zbrenty := b;
end;


{---------------------------------------------------------------------------}
function zbrent(f: TFuncX; a,b,t: extended; var ic,err: integer): extended;
  {-Brent/Dekker algorithm with guaranteed convergence for finding a zero  }
  { of a function: Return a zero x of the function f in the interval [a,b],}
  { function zbrenty with y=0.0 is used.                                   }
begin
  zbrent := zbrenty(f,0.0,a,b,t,ic,err);
end;


{---------------------------------------------------------------------------}
function zeroin(f: TFuncX; a,b,t: extended): extended;
  {-Return a zero x of the function f in the interval [a, b] to within a  }
  { tolerance 6*eps_x*|x| + 2*t, where t is a positive tolerance, assumes }
  { that f(a) and f(b) have different signs. Simplified version of zbrent.}
var
  ic,err: integer;
begin
  zeroin := zbrent(f,a,b,t,ic,err);
end;


{---------------------------------------------------------------------------}
function zridders(f: TFuncX; a,b,t: extended): extended;
  {-Find a zero x of f with Ridders' method in the interval [a,b] with}
  { a tolerance 2*eps_x*|x| + 0.5*t, returns NaN if f(a)*f(b) > 0.    }
var
  x1,x2,x3,xm,f1,f2,f3,fm,tol: extended;
  s1,s2,s3,sm: integer;
  icnt: integer;
const
  MAXIT = 100;
begin

  {Ref: Numerical recipes [13], 2nd ed, Ch. 9.2}
  zridders := Nan_x;

  x1 := a;
  f1 := f(x1);
  if f1=0.0 then begin
    zridders := x1;
    exit;
  end;

  x2 := b;
  f2 := f(x2);
  if f2=0.0 then begin
    zridders := x2;
    exit;
  end;

  s1 := isign(f1);
  s2 := isign(f2);
  s3 := s2-s1;

  {No bracket}
  if s3=0 then exit;

  x3   := Nan_x;
  icnt := 0;
  repeat
    inc(icnt);
    xm := (x1 + x2)/2;
    fm := f(xm);
    sm := isign(fm);
    if fm=0.0 then begin
      zridders := xm;
      exit;
    end;
    f3 := sqrt(sqr(fm) - f1 * f2);
    if f3<>0.0 then begin
      {See NR[13] 9.2.4}
      x3 := xm + (xm - x1) * isign(f1 - f2) * fm / f3;
      f3 := f(x3);
      s3 := isign(f3);
    end
    else begin
      zridders := x3;   {May be NaN, if very first denominator is zero}
      exit;
    end;
    if sm <> s3 then begin
      x1 := xm;  f1 := fm;  s1 := sm;
      x2 := x3;  f2 := f3;  s2 := s3;
    end
    else if s1 <> s3 then begin
      x2 := x3;  f2 := f3;  s2 := s3;
    end
    else if s2 <> s3 then begin
      x1 := x3;  f1 := f3;  s1 := s3;
    end
    else exit; {Should not happen, returns NaN}
    tol := 2.0*eps_x*abs(x3) + 0.5*t;
  until (abs(x1-x2) < tol) or (icnt > MAXIT);
  {$ifdef debug}
    if icnt > MAXIT then writeln('*** zridders: icnt > MAXIT');
  {$endif}
  {Return abscissa with smaller function value}
  if (abs(f1) < abs(f2)) then zridders := x1
  else zridders := x2;
end;



{---------------------------------------------------------------------------}
{---------------------- Quadratic/cubic equations --------------------------}
{---------------------------------------------------------------------------}

{Solve quadratic equations with double coefficients: Two functions   }
{implementing ideas of G.E. Forsythe, W. Kahan, P.H. Sterbenz (high  }
{precision calculation of discriminant, scaling by powers of two etc)}

{Formerly a separate unit, here the old version history:}

(*************************************************************************
 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     28.07.09  W.Ehrhardt  Initial BP7 version
 0.11     28.07.09  we          Other compilers ($N+}
 0.12     28.07.09  we          Selftest
 0.13     28.07.09  we          Fix x1<x2 if c=0
 0.14     18.08.09  we          squad_bsq
 0.15     18.08.09  we          fix squad_srs if c=0
 0.16     18.08.09  we          squad_s2 easy case, TP6+
 0.17     19.08.09  we          squad_s2 hard cases
 0.18     19.08.09  we          squad_s2 for TP5x
 0.19     19.08.09  we          improved TP5x ldexpd/getexpd
 0.20     19.08.09  we          improved squad_s2 case H3
 0.21     19.08.09  we          Handling of INFs and NANs
 0.22     07.09.09  we          Make y1 <= y2 in squad_core
 0.23     07.09.09  we          Return INFs in case H3
 0.24     08.09.09  we          squad_s2 renamed to squadx, code cleanup
 0.25     10.09.09  we          Bug fix ldexpd $ifdef VER5X
 0.26     27.01.10  we          uses amath, references renumbered
**************************************************************************)

{---------------------------------------------------------------------------}
function getexpd(d: double): longint;
  {-Return e with d=m*2^e and 1 <= abs(m) < 2}
var
  e: integer;
  x: extended;
begin
  if d=0 then getexpd := 0
  else begin
    x := d;
    e := THexExtW(x)[4] and $7FFF;
    if e=$7FFF then getexpd := e {INF}
    else getexpd := e-$3FFF;
  end;
end;


{---------------------------------------------------------------------------}
function IsInfOrNAN(a,b,c: double): boolean;
  {-Return true if at least one of a,b,c is +-INF or NAN}
begin
  IsInfOrNAN  := (THexDblW(a)[3] and $7FF0 = $7FF0) or
                 (THexDblW(b)[3] and $7FF0 = $7FF0) or
                 (THexDblW(c)[3] and $7FF0 = $7FF0);
end;


{---------------------------------------------------------------------------}
procedure break2(X: double; var Xh,Xt: double);
  {-Calculate Xh = X rounded to 26 sig. bits and Xt = X-Xh exactly in 26 sig.}
  { bits, so products like Xh*Xh, Xh*Xt, Xt*Xt can all be computed exactly.}
var
  bigX,difX: double;
begin
  bigX := X*134217729.0; { = X*(2^27 + 1)}
  difX := X - bigX;
  Xh   := difX + bigX;
  Xt   := X - Xh;
end;


{---------------------------------------------------------------------------}
function dscrmt(a,b,c: double): double;
  {-(high precision) calculation of the discriminant: dscrmt(a,b,c) = b*b-a*c}
var
  d,e: double;
  ah,at,bh,bt,ch,ct,p,q,dp,dq: double;
begin
  {http://www.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf}
  p := b*b;
  q := a*c;
  d := p-q;
  {d is now the (approximate) discriminant. The next two steps check}
  {whether d is adequately accurate. If not, recompute in high precision.}
  e := p+q;
  if 3.0*abs(d) < e  then begin
    {High precision calculation necessary. Split a,b,c}
    break2(a,ah,at);
    break2(b,bh,bt);
    break2(c,ch,ct);
    {DO NOT TRY TO COMBINE the following single steps. The correctness of the}
    {algorithm depends on the feature that every step is rounded to 53 bits.}
    {dp = ((bh*bh - p) + 2*bh*bt) + bt*bt}
    dp := bh*bh;
    dp := dp - p;
    d  := bh*bt;
    dp := dp + d;
    dp := dp + d;
    d  := bt*bt;
    dp := dp + d;
    {dq = ((ah*ch - q) + (ah*ct + at*ch)) + at*ct}
    dq := ah*ch;
    dq := dq - q;
    d  := ah*ct;
    dq := dq + d;
    d  := at*ch;
    dq := dq + d;
    d  := at*ct;
    dq := dq + d;
    d  := p - q;
    dp := dp - dq;
    d  := d  + dp;
  end;
  dscrmt := d;
end;


{---------------------------------------------------------------------------}
function squad_core(a,b,c: double; var x1,y1,x2,y2: double): integer;
  {-Solve the quadratic equation a*x^2 + b*x + c = 0. Result is number}
  { of different solutions: 0 (if a=b=0), 1 (x1), or 2 (x1,x2). If the}
  { result is = -2, x1+i*y1 and x2+i*y2 are the two complex solutions.}
  { No precautions against over/underflow, nor NAN/INF coefficients.}
var
  d: double;
begin

  {Setup default return values}
  squad_core := 2;
  x2 := 0.0;
  y1 := 0.0;
  y2 := 0.0;

  if a=0.0 then begin
    {solve bx+c=0}
    if b<>0.0 then begin
      squad_core := 1;
      x1 := -c/b;
    end
    else squad_core := 0;
    exit;
  end;

  if c=0.0 then begin
    {ax^2+bx = 0 = x(ax+b) with a<>0}
    x1 := -b/a;
    if x1>0.0 then begin
      x2 := x1;
      x1 := 0.0;
    end;
    exit;
  end;

  {Here a<>0 and c<>0. Divide b by -2 (without loss of precision) to simplify formulas}
  if b<>0.0 then b := -0.5*b;
  d := dscrmt(a,b,c);
  if d>0.0 then begin
    d := sqrt(d);
    if b >= 0.0 then d := b+d
    else d := b-d;
    x1 := d/a;
    x2 := c/d;
    {make x1 <= x2}
    if x1 > x2 then begin
      d  := x2;
      x2 := x1;
      x1 := d
    end;
  end
  else if d<0.0 then begin
    squad_core := -2;
    d  := sqrt(-d);
    x1 := b/a;
    x2 := x1;
    {make y1 <= y2}
    y2 := abs(d/a);
    y1 := -y2;
  end
  else begin
    {d = 0, real double root}
    squad_core := 1;
    x1 := b/a;
    x2 := x1;
  end;
end;


{---------------------------------------------------------------------------}
function squad(a,b,c: double; var x1,y1,x2,y2: double): integer;
  {-Solve the quadratic equation a*x^2 + b*x + c = 0. Result is the number}
  { of different solutions: 0 (if a=b=0), 1 (x1), or 2 (x1,x2). If the}
  { result is = -2, x1+i*y1 and x2+i*y2 are the two complex solutions.}
  { No precautions against over/underflow, NAN/INF coefficients return 0.}
begin
  if IsInfOrNAN(a,b,c) then squad := 0
  else squad := squad_core(a,b,c,x1,y1,x2,y2);
end;


{---------------------------------------------------------------------------}
function squadx(a,b,c: double; var x1,y1,x2,y2: double): integer;
  {-Solve the quadratic equation a*x^2 + b*x + c = 0. Result is the number}
  { of different solutions: 0 (if a=b=0 or INF/NAN), 1 (x1), or 2 (x1,x2).}
  { If the result is = -2, then x1 + i*y1 and x2 + i*y2 are the two complex}
  { solutions. Uses scaling by powers of two to minimize over/underflows.}
var
  ea,eb,ec,K,L,M: longint;
  t: double;
  sdiff: boolean;

const
  MaxExpDbl = 1019; {max. usable double exponent (for 4ac with 1<=|a|,|c|<2)}

  function setinf(x,y: double): double;
    {-Return sign(x)*sign(y)*INF}
  begin
    {return sign(x)*sign(y)*INF}
    if (THexDblW(x)[3] xor THexDblW(y)[3]) and $8000 = 0 then setinf := PosInf_d
    else setinf := NegInf_d;
  end;

begin
  if IsInfOrNAN(a,b,c) then begin
    squadx := 0;
    exit;
  end;

  if (a=0.0) or (c=0.0) then squadx := squad_core(a,b,c,x1,y1,x2,y2)
  else begin
    {The main idea is to scale x and the coefficients with separate scaling}
    {factors: Substitute x = 2^K*y and multiply the coefficients by 2^L.}
    {This is described (with base 16) in P. Sterbenz' book [17], which is out}
    {print since many years, an available paper based on [17] is W.Y. Sit [18]}

    {Get the coefficient exponents ea,eb,ec with a=2^ea*a', 1<=a'<2 etc}
    ea := getexpd(a);
    eb := getexpd(b);
    ec := getexpd(c);
    if b=0.0 then begin
      L  := -ea;
      K  := 0;
    end
    else begin
      L  := ea-eb-eb;
      K  := eb-ea;
    end;
    {c can only be scaled if new exponent ec+L is in double range}
    if abs(ec+L) <= MaxExpDbl then begin
      {Easy case, scale and perform standard procedure}
      a := ldexpd(a, -ea);
      b := ldexpd(b, -eb);
      c := ldexpd(c, L);
      squadx := squad_core(a,b,c,x1,y1,x2,y2);
      if K<>0 then begin
        {rescale results}
        x1 := ldexpd(x1, K);
        x2 := ldexpd(x2, K);
        y1 := ldexpd(y1, K);
        y2 := ldexpd(y2, K);
      end
    end
    else begin
      {Hard cases: c cannot be scaled. Remember if a*c < 0.}
      sdiff := ((a<0.0) and (c>0.0)) or ((a>0.0) and (c<0.0));
      {Setup default return values}
      squadx := 2;
      y1 := 0.0;
      y2 := 0.0;
      if b=0.0 then begin
        a := ldexpd(a, -ea);
        {here 1<=|a|<2, c cannot be scaled. Compute t = sqrt(abs(c/a)) }
        {scale c using an even exponent near exponent(sqrt(c))}
        M := (ec+L) div 2;
        c := ldexpd(c, L-M-M);
        t := ldexpd(sqrt(abs(c)/abs(a)),K+M);
        if sdiff then begin
          {Case H1: 4ac>0, roots are +- t, make x1 <= x2}
          x1 := -t;
          x2 := +t;
        end
        else begin
          {Case H2: 4ac>0, roots are +- t*i}
          {b^2 - 4ac < 0}
          x1 := 0.0;
          x2 := 0.0;
          y1 := -t;
          y2 := +t;
          squadx := -2;
        end;
      end
      else begin
        if ec+L < -MaxExpDbl then begin
          {Case H3: b<>0, |4ac| << b^2: ignore the 4ac term}
          {here direct overflows could occur depending on values of a,b,c}
          if eb-ea > MaxExpDbl then x1 := SetInf(-b,a) else x1 := -b/a;
          {x2 = (c/a)/x1 = -(c/a)/(b/a) = -(c/a)*(a/b) = -c/b}
          if ec-eb > MaxExpDbl then x2 := SetInf(-c,b) else x2 := -c/b;
          {make x1 <= x2}
          if x1 > x2 then begin
            t  := x2;
            x2 := x1;
            x1 := t;
          end;
        end
        else begin
          {here b<>0, ec+L > MaxExpDbl, |4ac| >> b^2}
          {Ignore b, calculate t as above}
          a := ldexpd(a, -ea);
          b := ldexpd(b, -eb);
          M := (ec+L) div 2;
          c := ldexpd(c, L-M-M);
          t := ldexpd(sqrt(abs(c)/abs(a)),K+M);
          if sdiff then begin
            {Case H4, real roots, make x1 <= x2}
            x1 := -t;
            x2 := +t;
          end
          else begin
            {Case H5: Complex roots with real part -b/(2a)}
            x1 := -0.5*b/a;
            if K<>0 then x1 := ldexpd(x1, K);
            x2 := x1;
            y1 := -t;
            y2 := +t;
            squadx := -2;
          end;
        end;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure cubeval(x,a,b,c,d: double; var q,q1,b1,c2: double);
  {-Internal use: return q = q(x) = ax^3 + bx^2 + cx + d,}
  { and q1 = q'(x), b1 = a*x + b, and  c2 := b1*x + c.}
var
  q0: double;
begin
  {Ref:  W.Kahan [60], procedure EVAL}
  q0 := a*x;
  b1 := q0 + b;
  c2 := b1*x + c;
  q1 := (q0+b1)*x + c2;
  q  := c2*x + d;
end;


{---------------------------------------------------------------------------}
procedure cubsolve(a,b,c,d: double; var x,x1,y1,x2,y2: double);
  {-Solve the cubic equation ax^3 + bx^2 + cx + d = 0: compute a real root}
  { x (may be INF if a~0) and two complex zeros x1 + i*y1, x2 + i*y2 where}
  { y2 = -y1 may be zero, i.e. there are three reel roots.}
var
  n: integer;
  q,q1,b1,c2,r,s,t,x0,s1: double;

const
  maxba = 0.225752123764894487e103;  {~cbrt(MaxDouble)/2.5}

begin
  if IsInfOrNAN(a,b,c) or IsNaNOrInf(d) then begin
    x  := NaN_d;
    x1 := NaN_d;
    x2 := NaN_d;
    y1 := 0.0;
    y2 := 0.0;
    exit;
  end;

  {Based on  W.Kahan [60], procedure QBC with two modifications:}
  {Handle very small a and sort 3 real roots in ascending order.}
  t  := abs(a);
  r  := abs(b);
  if r/maxba >= t then begin
    {a=0 or b/a very large , use x=-b/a and solve quadratic.}
    if (t=0.0) or ((t>=1.0) and (r >= MaxDouble/t)) then begin
      x := PosInf_d;
      n := -isign(a)*isign(b);
      if n<>0 then x := x*n;
    end
    else x := -b/a;
    {Avoid warning}
    if 0=squad(b,c,d,x1,y1,x2,y2) then ;
  end
  else if d=0.0 then begin
    x := 0.0;
    if 0=squad(a,b,c,x1,y1,x2,y2) then ;
  end
  else begin
    x := -(b/a)/3.0;
    cubeval(x,a,b,c,d,q,q1,b1,c2);
    t := q/a;
    r := cbrt(abs(t));
    s := isign(t);
    t := -q1/a;
    if t>0.0 then r := 1.324718*maxx(r, sqrt(t));
    x0 := x - s*r;
    if x<>x0 then begin
      s1 := succd(1.0);
      repeat
        x := x0;
        cubeval(x,a,b,c,d,q,q1,b1,c2);
        if q1=0.0 then x0 := x
        else x0 := x - (q/q1)/s1;
      until s*x0 <= s*x;
      if x<>0.0 then begin
        t := -d/x;
        if abs(a*x*x) > abs(t) then begin
          c2 := t;
          b1 := (c2-c)/x;
        end;
      end;
    end;
    if 0=squad(a,b1,c2,x1,y1,x2,y2) then ;
  end;
  {Sort 3 real roots in ascending order}
  if (y1=0.0) and (y2=0.0) then begin
    if x > x1 then begin
      {swap x,x1}
      t  := x;
      x  := x1;
      x1 := t;
    end;
    if x1 > x2 then begin
      {swap x1,x2}
      t  := x1;
      x1 := x2;
      x2 := t;
      if x > x1 then begin
        {swap x,x1}
        t  := x;
        x  := x1;
        x1 := t;
      end;
    end;
  end;
end;


type
  TRVMatrix  = array[1..MaxDeg, 1..MaxDeg] of double;

{---------------------------------------------------------------------------}
procedure balanc(var a: TRVMatrix; n: integer);
  {-Balance n x n submatrix of a, based on EISPACK balanc }
  { This is special version for the use of hrq computing  }
  { roots of polynomials. No scale information is returned}
  { and almost always I found LOW=1 and IGH=1. It is a    }
  { port from [13], balanc.c}
const
  radix=2.0;
var
  last,j,i: integer;
  s,r,g,f,c,sqrdx: double;
begin
  sqrdx := sqr(radix);
  last := 0;
  while last=0 do begin
    last := 1;
    {Calculate row and column norms}
    for i:=1 to n do begin
      c := 0.0;
      r := 0.0;
      for j:=1 to n do begin
        if j<>i then begin
          c := c + abs(a[j,i]);
          r := r + abs(a[i,j])
        end;
      end;
      {guard against zero c or r due to underflow}
      if (c<>0.0) and (r<>0.0) then begin
        {If both are nonzero,}
        g := r/radix;
        f := 1.0;
        s := c + r;
        {find the integer power of the machine radix}
        {that comes closest to balancing the matrix }
        while c<g do begin
          f := f*radix;
          c := c*sqrdx
        end;
        g := r*radix;
        while c>g do begin
          f := f/radix;
          c := c/sqrdx
        end;
        {now balance}
        if (c+r)/f < 0.95*s then begin
          last := 0;
          g := 1.0/f;
          {Apply similarity transformation}
          for j:=1 to n do a[i,j] := a[i,j]*g;
          for j:=1 to n do a[j,i] := a[j,i]*f
        end
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure hqr(
  n: integer;          {actual dimension of the square matrix  }
  blow: integer;       {org low, from general balance or 1     }
  bhigh: integer;      {org igh, from general balance or n     }
  var h: TRVMatrix ;   {the Matrix, the n x n part is processed}
  var wr,wi: TPolyVec; {real, imaginary eigenvalues.           }
  var ierr: integer);  {error code, see below                  }

  {-Find the eigenvalues of a real upper Hessenberg matrix by the QR method}

var
  en,enm2,i,itn,its,j,k,l,ll,m,mm,mp2,na: integer;
  notlas: boolean;
  p,q,r,s,t,tst1,tst2,w,x,y,zz,norm: double;
label
  60,70,100,150,260,270,280;

{---------------------------------------------------------------------------}
{----------  Original comment from the EISPACK routime hqr.f  --------------}
{---------------------------------------------------------------------------}
(*
c  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
c
c
c     this subroutine is a translation of the algol procedure hqr,
c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
c
c     this subroutine finds the eigenvalues of a real
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        h contains the upper hessenberg matrix.  information about
c          the transformations used in the reduction to hessenberg
c          form by  elmhes  or  orthes, if performed, is stored
c          in the remaining triangle under the hessenberg matrix.
c
c     on output
c
c        h has been destroyed.  therefore, it must be saved
c          before calling  hqr  if subsequent calculation and
c          back transformation of eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated september 1989.
c
c     ------------------------------------------------------------------
*)

begin
     {Ref: http://www.netlib.org/eispack/hqr.f}
     {The next statements makes some compilers happy}
     p := 0.0;
     q := 0.0;
     r := 0.0;
     l := 0;
     m := 0;

     {start original code}
     ierr := 0;
     norm := 0.0;
     k := 1;

     {store roots isolated by balanc and compute matrix norm}
     for i:=1 to n do begin
       for j:=k to n do norm := norm+abs(h[i,j]);
       k := i;
       if (i < blow) or (i > bhigh) then begin
         wr[i] := h[i,i];
         wi[i] := 0.0;;
       end;
     end;

     en := bhigh;
     t := 0.0;
     itn := 30*n;

     {search for next eigenvalues}
60:  if en<blow then exit;
     its := 0;
     na := en-1;
     enm2 := na-1;

     {look for single small sub-diagonal element}
     {         for l=en step -1 until low do -- }
     {*WE: l must be valid if loop finished     }
70:  for ll:=blow to en do begin
       l := en + blow - ll;
       if l=blow then goto 100;
       s := abs(h[l-1,l-1]) + abs(h[l,l]);
       if s=0.0 then s := norm;
       tst1 := s;
       tst2 := tst1 + abs(h[l,l-1]);
       if tst2=tst1 then goto 100;
     end;

     {form shift}
100: x := h[en,en];
     if l=en then goto 270;
     y := h[na,na];
     w := h[en,na]*h[na,en];
     if l=na then goto 280;

     if itn=0 then begin
       ierr := en;
       exit;
     end;

     if (its=10) or (its=20) then begin
       {form exceptional shift}
       t := t+x;
       for i:=blow to en do h[i,i] := h[i,i] - x;
       s := abs(h[en,na]) + abs(h[na,enm2]);
       x := 0.75*s;
       y := x;
       w := -0.4375*s*s;
     end;

     its := its+1;
     itn := itn-1;

     {look for two consecutive small sub-diagonal elements.}
     {            for m=en-2 step -1 until l do --}
     {*WE: m must be valid if loop finished}
     for mm:=l to enm2 do begin
       m := enm2 + l - mm;
       zz:= h[m,m];
       r := x - zz;
       s := y - zz;
       p := (r*s - w)/h[m+1,m] + h[m,m+1];
       q := h[m+1,m+1] - zz - r - s;
       r := h[m+2,m+1];
       s := abs(p) + abs(q) + abs(r);
       p := p/s;
       q := q/s;
       r := r/s;
       if m=l then goto 150;
       tst1 := abs(p)*(abs(h[m-1,m-1]) + abs(zz) +abs(h[m+1,m+1]));
       tst2 := tst1 + abs(h[m,m-1])*(abs(q) + abs(r));
       if tst2=tst1 then goto 150;
     end;

150: mp2 := m+2;

     for i:=mp2 to en do begin
       h[i,i-2] := 0.0;
       if i<>mp2 then h[i,i-3] := 0.0;
     end;

     {double qr step involving rows l to en and columns m to en}
     for k:=m to na do begin
       notlas := k <> na;
       if k<>m then begin
         p := h[k,k-1];
         q := h[k+1,k-1];
         r := 0.0;
         if notlas then r := h[k+2,k-1];
         x := abs(p) + abs(q) + abs(r);
         if x=0.0 then goto 260;
         p := p/x;
         q := q/x;
         r := r/x;
       end;

       s := copysignd(sqrt(p*p + q*q + r*r),p);
       if k<>m then h[k,k-1] := -s*x
       else if l<>m then h[k,k-1] := -h[k,k-1];

       p := p+s;
       x := p/s;
       y := q/s;
       zz:= r/s;
       q := q/p;
       r := r/p;

       if not notlas then begin
         {row modification}
         for j:=k to en do begin
           p := h[k,j] + q*h[k+1,j];
           h[k,j] := h[k,j] - p*x;
           h[k+1,j] := h[k+1,j] - p*y;
         end;

         {j = min0(en,k+3)}
         j := k+3;
         if en<j then j:= en;
         {column modification}
         for i:=l to j do begin
           p := x*h[i,k] + y*h[i,k+1];
           h[i,k] := h[i,k] - p;
           h[i,k+1] := h[i,k+1] - p*q;
         end;
       end
       else begin
         {row modification}
         for j:=k to en do begin
           p := h[k,j] + q*h[k+1,j] + r*h[k+2,j];
           h[k,j] := h[k,j]-p*x;
           h[k+1,j] := h[k+1,j] - p*y;
           h[k+2,j] := h[k+2,j] - p*zz;
         end;
         {column modification}
         {j = min0(en,k+3)}
         j := k+3;
         if en<j then j:= en;
         for i:=l to j do begin
           p := x*h[i,k] + y*h[i,k+1] + zz*h[i,k+2];
           h[i,k] := h[i,k]-p;
           h[i,k+1] := h[i,k+1] - p*q;
           h[i,k+2] := h[i,k+2] - p*r;
         end;
       end;
260: end;

     goto 70;

     {one root found}
270: wr[en] := x + t;
     wi[en] := 0.0;
     en := na;
     goto 60;

     {two roots found}
280: p := (y-x)/2.0;
     q := p*p + w;
     zz := sqrt(abs(q));
     x := x+t;
     if q >=0.0 then begin
       {real pair}
       zz := p + copysignd(zz,p);
       wr[na] := x+zz;
       wr[en] := wr[na];
       if zz <> 0.0 then wr[en] := x-w/zz;
       wi[na] := 0.0;
       wi[en] := 0.0;
     end
     else begin
       {complex pair}
       wr[na] := x+p;
       wr[en] := x+p;
       {*WE: put root with negative imaginary first}
       wi[na] := -zz;
       wi[en] := zz;
     end;

     en := enm2;
     goto 60;
end;


{---------------------------------------------------------------------------}
procedure PolyRootSort(n: integer; var x,y: TPolyVec);
  {-Sort the roots x[k] + i*y[k] with insertion-sort,}
  { assumes that a complex pair is already sorted.   }
var
  i,j: integer;
  d: double;
begin
  for i:=2 to n do begin
    j := i;
    while (j>1) and (x[j-1] > x[j]) do begin
      d := x[j];
      x[j] := x[j-1];
      x[j-1] := d;
      d := y[j];
      y[j] := y[j-1];
      y[j-1] := d;
      j := j-1;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure PolyRootPolish(const p: TPolyVec; n: integer; var x: double);
  {-Do one Newton step for p(x)}
var
  px,dp,z: extended;
begin
  {compute px = p(x) and dp = p'(x)}
  PolEvalDeriv(x,p,n+1,px,dp);
  if (abs(px) < 2*eps_d) or (abs(dp) < 8*eps_d) then exit;
  {step}
  z := x - px/dp;
  if abs(PolEval(z,p,n+1)) < abs(px) then x := z;
end;


{---------------------------------------------------------------------------}
procedure PolyRootPolishComplex(const p,dp: TPolyVec; n: integer; var x,y: double);
  {-Do one Newton step for p(x), dp = p('x)}
var
  x0,x1,y0,y1: double;
  dx,dy,a0: double;
begin
  {exit if p(x+iy) = 0}
  PolEvalC(p,n+1,x,y,x0,y0);
  a0 := hypot(x0,y0);
  if a0=0.0 then exit;

  {exit if p'(x+iy) ~ 0}
  PolEvalC(dp,n,x,y,x1,y1);
  if hypot(x1,y1) < 2*eps_d then exit;

  {complex division p(x+iy)/p'(x+iy)}
  dx := (x0*x1 + y0*y1)/(x1*x1 + y1*y1);
  dy := (x1*y0 - y1*x0)/(x1*x1 + y1*y1);

  {Compute new x,y}
  x1 := x - dx;
  y1 := y - dy;

  PolEvalC(p,n+1,x1,y1,x0,y0);
  if hypot(x0,y0) < a0 then begin
    {update only if new |p(x+iy)| is smaller}
    x := x1;
    y := y1;
  end;
end;


{---------------------------------------------------------------------------}
procedure PolyRoots(const p: TPolyVec; n: integer; var x,y: TPolyVec; var ierr: integer);
  {-Compute the n (complex) roots x[k] + i*y[k] of the polynomial  }
  { p(z) = p[0] + p[1]*z + ... p[n]*z^n, n>0, p[n]<>0. Real roots  }
  { have y[k]=0. ierr is an error code 0: OK, -1: n<1, -2: p[n]=0, }
  { k: iteration count exceeded while the kth root is being sought.}
  { n must be <= MaxDeg (currently = 20); the cases n=1 and 2 are  }
  { handled separately: n=1 directly and n=2 with squad,           }

  { PolyRoots uses a companion matrix method, balancing, and the   }
  { QR algorithm for the eigenvalues of an upper Hessenberg matrix.}
var
  ph: ^TRVMatrix;
  a: TPolyVec;
  i,m: integer;
begin
  if (n<1) or (n>MAXDeg) then begin
    ierr := -1;
    exit;
  end;
  if p[n]=0.0 then begin
    ierr := -2;
    exit;
  end;

  {zero-fill the roots}
  for i:=0 to n do begin
    x[i] := 0;
    y[i] := 0;
  end;

  ierr := 0;

  {detect low power zero coefficients}
  if p[0]=0.0 then begin
    i := 1;
    while p[i]=0.0 do inc(i);
    m := -1;
    while i<=n do begin
      inc(m);
      a[m] := p[i];
      inc(i);
    end;
  end
  else begin
    m := n;
    a := p;
  end;

  if m<3 then begin
    {case degree < 3}
    if m=1 then begin
      {a[0] + a[1]*z = 0}
      x[1] := -a[0]/a[1];
      y[1] := 0.0;
    end
    else begin
      {error -3 if squad cannot find any root}
      if squad(a[2], a[1], a[0], x[1], y[1], x[2], y[2])=0 then ierr := -3;
    end;
  end
  else begin
    {case degree > 2}
    {allocate and zerofill the matrix}
    ph := calloc(sizeof(TRVMatrix));
    if ph=nil then begin
      ierr := -4;
      exit;
    end;
    {setup the companion matrix}
    for i:=1 to m do ph^[1,i]   := -a[m-i]/a[m];
    for i:=2 to m do ph^[i,i-1] := 1.0;
    {balance the matrix}
    Balanc(ph^,m);
    {compute eigenvalues/roots}
    HQR(n,1,m,ph^,x,y,ierr);
    {release memory}
    cfree(ph, sizeof(TRVMatrix));
  end;
end;


{---------------------------------------------------------------------------}
procedure PolyRootsA(const p: TPolyVec; n: integer; var x,y: TPolyVec; var ierr: integer);
  {-Perform PolyRoots, then improve and sort the roots}
var
  i: integer;
  dp: TPolyVec; {p'(x)}
begin
  PolyRoots(p,n,x,y,ierr);
  if ierr=0 then begin
    {build p'(x)}
    dp[n]:=0;
    for i:=n-1 downto 0 do dp[i] := (i+1)*p[i+1];
    {Polish non-zero roots}
    for i:=1 to n do begin
      if (y[i]<>0.0) or (x[i]<>0.0) then begin
        if y[i]=0.0 then begin
          {Polish real root}
          PolyRootPolish(p, n, x[i])
        end
        else begin
          {Polish complex root}
          PolyRootPolishComplex(p, dp, n, x[i], y[i]);
        end;
      end;
    end;
    {sort roots, including trivial zeroes}
    PolyRootSort(n,x,y);
  end;
end;


{---------------------------------------------------------------------------}
procedure PolyRootsOA(const p: array of double; n: integer; var x,y: TPolyVec; var ierr: integer);
  {-Open array version, n <= high(p): perform PolyRoots, then improve and sort the roots}
var
  k: integer;
  coef: TPolyVec;
begin
  if (n > high(p)) or (n > MaxDeg) then begin
    ierr := -1;
    exit;
  end;
{$ifdef FPC}
  {Suppress warnings: Local variable does not seem to be initialized}
  coef[0] := 0.0;
{$endif}
  for k:=0 to n do coef[k] := p[k];
  PolyRootsA(coef,n,x,y,ierr);
end;


{---------------------------------------------------------------------------}
{----------------------- Convergence acceleration --------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure levinu1(an: extended; n,nmax: integer;
                  var qnum, qden: array of extended;
                  var extsum, sn: extended; var ierr: integer);
  {-Perform one step of the Levin u-transformation for sum(an,n=0..). }
  { The driver routine has to analyze the convergence of the process. }
  { DO NOT CHANGE the var parameters between successive calls!        }
  {#F}
  { an     - the n'th term of the sum                                 }
  { n      - the index of current term for this step                  }
  { nmax   - maximum valid value for n and index the a arrays         }
  { qnum   - backward diagonal of numerator array,   at least 0..nmax }
  { qden   - backward diagonal of denominator array, at least 0..nmax }
  { extsum - extrapolated value of the sum                            }
  { sn     - the n'th partial sum, will be set to an if n=0           }
  { ierr   - error if <> 0 (currently for n < 0,  n > nmax)           }
  {#F}
var
  j,k: integer;
  t,f1,factor,ratio: extended;
begin

  {Ref: Fessler et.al. [57] and Alg602, subroutine WHIZ1}
  {See also: Weniger [56], Ch.7 and subroutine GLEVIN   }
  ierr := 0;
  if (n<0) or (n>nmax) then begin
    ierr := -1;
    exit;
  end;

  if n=0 then sn := an
  else sn := an + sn;

  k  := n;
  f1 := k+1.0;
  t  := an*sqr(f1);
  if abs(t)<Sqrt_MinExt then begin
    qden[k] := 1.0;
    qnum[k] := 1.0;
  end
  else begin
    qden[k] := 1.0/t;
    qnum[k] := sn*qden[k];
  end;
  if k>0 then begin
    factor := 1.0;
    ratio  := k/f1;
    for j:=k-1 downto 0 do begin
      t := factor*(j+1)/f1;
      factor  := factor*ratio;
      qden[j] := qden[j+1] - t*qden[j];
      qnum[j] := qnum[j+1] - t*qnum[j];
    end;
  end;
  {If quotient will overflow, return large value}
  t := abs(qden[0]);
  if (t>1.0) or (t*MaxExtended > abs(qnum[0])) then begin
    extsum := qnum[0]/qden[0];
  end
  else begin
    extsum := Sqrt_MaxExt;
  end;
end;


{---------------------------------------------------------------------------}
procedure wynneps1(an: extended; n,nmax: integer; sum: boolean;
                   var e: array of extended;
                   var extlim, sn: extended; var ierr: integer);
  {-Perform one step of Wynn's epsilon algorithm for the sequence  }
  { an, n=0.. or the sum(an,n=0..). Note: The driver routine has to}
  { analyze the convergence of the whole process. DO NOT CHANGE the}
  { var parameters between successive calls!                       }
  {#F}
  { an     - the n'th term ot the sequence or sum                  }
  { n      - the index of current term for this step               }
  { nmax   - maximum valid value for n and index the e array       }
  { sum    - if true then a sum is processed, otherwise a sequence }
  { e      - working array for epsilon table, at least 0..nmax     }
  { extlim - extrapolated value of the sequence or sum             }
  { sn     - the n'th partial sum, will be set to an if n=0        }
  { ierr   - error if <> 0 (currently for n < 0,  n > nmax)        }
  {#F}
var
  aux1,aux2,diff: extended;
  j: integer;
begin

  {Ref: Weniger [56], Ch.4 and subroutine EPSAL}
  ierr := 0;
  if (n<0) or (n>nmax) then begin
    ierr := -1;
    exit;
  end;

  if sum and (n>0) then sn := sn + an
  else sn := an;

  e[n] := sn;
  if n=0 then extlim := sn
  else begin
    aux2 := 0.0;
    for j:=n downto 1 do begin
       aux1 := aux2;
       aux2 := e[j-1];
       diff := e[j] - aux2;
       if abs(diff) <= Sqrt_MinExt then e[j-1] := aux1
       else e[j-1] := aux1 + 1.0/diff;
    end;
    extlim := e[n and 1];
  end;
end;



{---------------------------------------------------------------------------}
{--------------------- Numerical integration -------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure quanc8(fun: TFuncX; a,b,abserr,relerr: extended;
                 var result, errest, flag: extended; var nofun: longint);
  {-Estimate the integral of fun(x) from a to b to a user provided tolerance.}
  { Pascal translation of the Fortran subroutine by Forsythe, Malcolm, Moler.}

  {An automatic adaptive routine based on the 8-panel Newton-Cotes rule: }
  {Input:                                                                }
  {  fun     The name of the integrand function subprogram fun(x).       }
  {  a       The lower limit of integration.                             }
  {  b       The upper limit of integration, b may be less than a.       }
  {  relerr  A  relative error tolerance,  should be non-negative.       }
  {  abserr  An absolute error tolerance,  should be non-negative.       }
  {Output:                                                               }
  {  result  An approximation to the integral hopefully satisfying the   }
  {          least stringent of the two error tolerances.                }
  {  errest  An estimate of the magnitude of the actual error.           }
  {  nofun   The number of function values used in calculation of result.}
  {  flag    A reliability indicator: if flag is zero, then result       }
  {          probably satisfies the error tolerance. If flag is          }
  {          xxx.yyy , then xxx = the number of intervals which have.    }
  {          not converged and 0.yyy = the fraction of the interval      }
  {          left to do when the limit on nofun was approached.          }
  {                                                                      }
  {Quanc8 is based on piecewise polynomial approximation and hence is not}
  {designed to handle certain kinds of integrals. Roughly, these are     }
  {integrals of functions for which some derivative up to 10th order is  }
  {unbounded or fails to exist.                                          }

var
  area, stone, step, cor11, temp: extended;
  qprev, qnow, qdiff, qleft, esterr, tolerr: extended;
  levmax, lev, i, j: integer;
  nim, nofin: longint;
  qright: array[1..31] of extended;
  f, x: array[0..16] of extended;
  fsave, xsave: array[1..8,1..30] of extended;

label
  stage3, stage5, stage6, stage7, stage8;

const
  levmin = 1;
  levout = 6;
  nomax  = 5000;

const
  w0: extended =   3956 / 14175; {Newton-Cotes weights multiplied by 8}   {Fix311}
  w1: extended =  23552 / 14175;
  w2: extended =  -3712 / 14175;
  w3: extended =  41984 / 14175;
  w4: extended = -18160 / 14175;

begin

  {This is my Pascal translation of subroutine QUANC8 from [29], Ch. 5.5.}
  {The overall structure of the Fortran code is left unchanged including }
  {the stage labels. Since the Pascal code should run on versions without}
  {'break' and 'continue', structured programming could be achieved only }
  {with artificial booleans or case variables.}

  {*** stage 1 *** general initialization}
  {set 'constants'}
  levmax := 30;

  {trouble when nofun reaches nofin}
  nofin  := nomax - 8*(levmax - levout + (2 shl levout));

  {initialize running sums to zero}
  flag   := 0.0;
  result := 0.0;
  cor11  := 0.0;
  errest := 0.0;
  area   := 0.0;
  nofun  := 0;
  if a=b then exit;

{$ifdef FPC}
  {Suppress warnings: Local variable does not seem to be initialized}
  fsave[1,1] := 0.0;
  xsave[1,1] := 0.0;
{$endif}

  {WE modification: make tolerances non-negative}
  if abserr < 0.0 then abserr := 0.0;
  if relerr < 0.0 then relerr := 0.0;

  {*** stage 2 *** initialization for first interval}
  lev   := 0;
  nim   := 1;
  qprev := 0.0;
  stone := (b - a) / 16.0;

  x[0]  := a;
  x[16] := b;
  x[8 ] := 0.5*(x[ 0] + x[16]);
  x[4 ] := 0.5*(x[ 0] + x[ 8]);
  x[12] := 0.5*(x[ 8] + x[16]);
  x[2 ] := 0.5*(x[ 0] + x[ 4]);
  x[6 ] := 0.5*(x[ 4] + x[ 8]);
  x[10] := 0.5*(x[ 8] + x[12]);
  x[14] := 0.5*(x[12] + x[16]);

  for i:=0 to 8 do begin
    j    := 2*i;
    f[j] := fun(x[j]);
    inc(nofun);
  end;

stage3: {*** central calculation ****}
  {requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16}
  {calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area}
  for i:=0 to 7 do begin
    j    := 2*i+1;
    x[j] := 0.5*(x[j-1] + x[j+1]);
    f[j] := fun(x[j]);
    inc(nofun);
  end;
  step  := (x[16] - x[0])/16.0;
  qleft := (w0*(f[0]+f[8]) + w1*(f[1]+f[7]) + w2*(f[2]+f[6]) +
          + w3*(f[3]+f[5]) + w4*f[4])*step;
  qright[lev+1] := (w0*(f[8]+f[16]) + w1*(f[9]+f[15]) + w2*(f[10]+f[14])
                  + w3*(f[11]+f[13]) + w4*f[12])*step;
  qnow  := qleft + qright[lev+1];
  qdiff := qnow  - qprev;
  area  := area  + qdiff;

  {*** stage 4 *** interval convergence test}
  esterr := abs(qdiff) / 1023.0;
  tolerr := maxx(abserr,relerr*abs(area)) * (step/stone);

  if lev <  levmin then goto stage5;
  if lev >= levmax then begin
    {current level is levmax}
    flag := flag + 1.0;
    goto stage7;
  end;
  if nofun  >  nofin  then goto stage6;
  if esterr <= tolerr then goto stage7;

stage5: {*** no convergence ***}
  {locate next interval}
  nim := 2*nim;
  lev := lev+1;
  {store right hand elements for future use}
  for i:=1 to 8 do begin
    fsave[i,lev] := f[i+8];
    xsave[i,lev] := x[i+8];
  end;
  {assemble left hand elements for immediate use}
  qprev := qleft;
  for i:=8 downto 1 do begin
    j    := 2*i;
    f[j] := f[i];
    x[j] := x[i];
  end;
  goto stage3;

stage6: {*** trouble section ***}
  {number of function values is about to exceed limit}
  nofin  := 2*nofin;
  levmax := levout;
  flag   := flag + (b - x[0]) / (b - a);

stage7: {*** interval converged ***}
  {add contributions into running sums.}
  result := result + qnow;
  errest := errest + esterr;
  cor11  := cor11  + qdiff / 1023.0;
  {locate next interval}
  while odd(nim) do begin
    nim := nim div 2;
    dec(lev);
  end;
  inc(nim);
  if lev <= 0 then goto stage8;

  {assemble elements required for the next interval}
  qprev := qright[lev];
  x[0]  := x[16];
  f[0]  := f[16];
  for i:=1 to 8 do begin
    f[2*i] := fsave[i,lev];
    x[2*i] := xsave[i,lev];
  end;
  goto stage3;

stage8: {*** finalize and return ***}
  result := result + cor11;
  {make sure errest not less than roundoff level}
  if errest > 0 then begin
    temp := abs(result);
    while temp+errest=temp do errest := 2.0*errest;
  end;
end;


{---------------------------------------------------------------------------}
procedure qk21(f: TFuncX; a,b: extended; var result,abserr,resabs,resasc: extended);
  {-Integrate function f over (a,b) with 21-point Gauss-Kronrod rule}
  { a      - lower limit of integration                     }
  { b      - upper limit of integration                     }
  { result - approximation to the integral I of f over (a,b)}
  { abserr - estimate of the modulus of the absolute error  }
  { resabs - approximation to the integral |f| over (a,b)   }
  { resasc - approximation integral |f-I/(b-a)| over (a,b)  }
const
  wg: array[1..5] of extended = (
         0.066671344308688137593568809893332,  {weights of the 10-point gauss rule}
         0.149451349150580593145776339657697,
         0.219086362515982043995534934228163,
         0.269266719309996355091226921569469,
         0.295524224714752870173892994651338);
  xgk: array[1..11] of extended = (            {abscissae of the 21-point Gauss/Kronrod rule}
         0.995657163025808080735527280689003,
         0.973906528517171720077964012084452,
         0.930157491355708226001207180059508,
         0.865063366688984510732096688423493,
         0.780817726586416897063717578345042,
         0.679409568299024406234327365114874,
         0.562757134668604683339000099272694,
         0.433395394129247190799265943165784,
         0.294392862701460198131126603103866,
         0.148874338981631210884826001129720,
         0.0);
  wgk: array[1..11] of extended = (            {weights of the 10-point gauss rule}
         0.011694638867371874278064396062192,
         0.032558162307964727478818972459390,
         0.054755896574351996031381300244580,
         0.075039674810919952767043140916190,
         0.093125454583697605535065465083366,
         0.109387158802297641899210590325805,
         0.123491976262065851077958109831074,
         0.134709217311473325928054001771707,
         0.142775938577060080797094273138717,
         0.147739104901338491374841515972068,
         0.149445554002916905664936468389821);
var
  hlgth: extended; {half-length of the interval}
  centr: extended; {mid point of the interval}
  dabsc: extended; {abscissa}
  resg : extended; {result of the 10-point Gauss formula}
  resk : extended; {result of the 21-point Kronrod formula}
  reskh: extended; {approximation to the mean value of f over (a,b)}
  fval1,
  fval2: extended; {function values}
  fc,fsum,dhlgth, t: extended;
  fv1,fv2: array[1..10] of extended;
  j,k: integer;
begin
  {Ref: Quadpack[34], subroutine dqk21 in dqk21.f}
  centr  := 0.5*(b+a);
  hlgth  := 0.5*(b-a);
  dhlgth := abs(hlgth);
  {compute the 21-point Kronrod approximation to}
  {the integral, and estimate the absolute error}
  resg   := 0.0;
  fc     := f(centr);
  resk   := wgk[11]*fc;
  resabs := abs(resk);
  for j:=1 to 5 do begin
    k := j*2;
    dabsc  := hlgth*xgk[k];
    fval1  := f(centr-dabsc);
    fval2  := f(centr+dabsc);
    fv1[k] := fval1;
    fv2[k] := fval2;
    fsum   := fval1  + fval2;
    resg   := resg   + wg[j]*fsum;
    resk   := resk   + wgk[k]*fsum;
    resabs := resabs + wgk[k]*(abs(fval1)+abs(fval2));
  end;
  for j:=1 to 5 do begin
    k := j*2-1;
    dabsc  := hlgth*xgk[k];
    fval1  := f(centr-dabsc);
    fval2  := f(centr+dabsc);
    fv1[k] := fval1;
    fv2[k] := fval2;
    fsum   := fval1  + fval2;
    resk   := resk   + wgk[k]*fsum;
    resabs := resabs + wgk[k]*(abs(fval1)+abs(fval2));
  end;
  reskh  := 0.5*resk;
  resasc := wgk[11]*abs(fc-reskh);

  for j:=1 to 10 do resasc := resasc + wgk[j]*(abs(fv1[j]-reskh)+abs(fv2[j]-reskh));

  result := resk*hlgth;
  resabs := resabs*dhlgth;
  resasc := resasc*dhlgth;
  abserr := abs((resk-resg)*hlgth);
  {scale error}
  if (resasc<>0.0) and (abserr<>0.0) then begin
    t := sqrt(abs(abserr/resasc));
    abserr := resasc*minx(1.0,200.0*t*t*t)
  end;
  t := 50.0*eps_x;
  if resabs > MinExtended/t then abserr := maxx(t*resabs,abserr);
end;


{---------------------------------------------------------------------------}
procedure qk15i(f: TFuncX; bound: extended; inf: integer; a, b: extended;
                var result, abserr, resabs, resasc: extended);
  {-Integration with 15-point transformed Gauss-Kronrod rule}
const
  wg : array[1..8] of extended = ( {weights of the 7-point Gauss rule}
         0.0, 0.129484966168869693270611432679082,
         0.0, 0.279705391489276667901467771423780,
         0.0, 0.381830050505118944950369775488975,
         0.0, 0.417959183673469387755102040816327);
  xgk: array[1..8] of extended = (   {abscissae of the 15-point Kronrod rule}
         0.991455371120812639206854697526329,
         0.949107912342758524526189684047851,
         0.864864423359769072789712788640926,
         0.741531185599394439863864773280788,
         0.586087235467691130294144838258730,
         0.405845151377397166906606412076961,
         0.207784955007898467600689403773245,
         0.0);
  wgk: array[1..8] of extended = ( {weights of the 15-point Kronrod rule}
         0.022935322010529224963732008058970,
         0.063092092629978553290700663189204,
         0.104790010322250183839876322541518,
         0.140653259715525918745189590510238,
         0.169004726639267902826583426598550,
         0.190350578064785409913256402421014,
         0.204432940075298892414161999234649,
         0.209482141084727828012999174891714);
var
  centr: extended;  {mid point of the interval}
  absc,
  absc1,
  absc2: extended;  {abscissa}
  fval1,
  fval2: extended;  {function value}
  hlgth: extended;  {half-length of the interval}
  resg:  extended;  {result of the 7-point gauss formula}
  resk:  extended;  {result of the 15-point Kronrod formula}
  reskh: extended;  {approximation to the mean value of the transformed integrand}
  tabsc1,
  tabsc2: extended; {transformed abscissa}
  dinf, fc, fsum, t: extended;
  fv1, fv2: array[1..7] of extended;
  j: integer;
begin
  {Ref: Quadpack[34], subroutine dqk15i in dqk15i.f}
  dinf   := minx(1.0,inf);
  centr  := 0.5*(a+b);
  hlgth  := 0.5*(b-a);
  tabsc1 := bound + dinf*(1.0-centr)/centr;
  fval1  := f(tabsc1);

  if inf=2 then fval1 := fval1+f(-tabsc1);
  fc := (fval1/centr)/centr;

  {compute the 15-point Kronrod approximation to}
  {the integral, and estimate the error}
  resg   := wg[8]*fc;
  resk   := wgk[8]*fc;
  resabs := abs(resk);

  for j:=1 to 7 do begin
    absc   := hlgth*xgk[j];
    absc1  := centr - absc;
    absc2  := centr + absc;
    tabsc1 := bound + dinf*(1.0-absc1)/absc1;
    tabsc2 := bound + dinf*(1.0-absc2)/absc2;
    fval1  := f(tabsc1);
    fval2  := f(tabsc2);
    if inf=2 then begin
      fval1 := fval1+f(-tabsc1);
      fval2 := fval2+f(-tabsc2);
    end;
    fval1  := (fval1/absc1)/absc1;
    fval2  := (fval2/absc2)/absc2;
    fv1[j] := fval1;
    fv2[j] := fval2;
    fsum   := fval1  + fval2;
    resg   := resg   + wg[j]*fsum;
    resk   := resk   + wgk[j]*fsum;
    resabs := resabs + wgk[j]*(abs(fval1)+abs(fval2))
  end;
  reskh  := resk*0.5;
  resasc := wgk[8]*abs(fc-reskh);

  for j:=1 to 7 do resasc := resasc + wgk[j]*(abs(fv1[j]-reskh)+abs(fv2[j]-reskh));

  result := resk*hlgth;
  resasc := resasc*hlgth;
  resabs := resabs*hlgth;
  abserr := abs((resk-resg)*hlgth);

  {scale error}
  if (resasc<>0.0) and (abserr<>0.0) then begin
    t := sqrt(abs(abserr/resasc));
    abserr := resasc*minx(1.0,200.0*t*t*t)
  end;
  t := 50.0*eps_x;
  if resabs > MinExtended/t then abserr := maxx(t*resabs,abserr);
end;


{---------------------------------------------------------------------------}
procedure qpsrt(limit: integer; var last, maxerr: integer; var ermax: extended;
                var elist: TQXLimArr; var iord: TQILimArr; var nrmax: integer);
  {-Maintain the descending ordering in the list of the local error estimates}

  { limit  - maximum number of error estimates the list can contain}
  { last   - number of error estimates currently in the list}
  { maxerr - maxerr points to the nrmax-th largest error estimate currently in the list}
  { ermax  - nrmax-th largest error estimate ermax = elist[maxerr]}
  { elist  - vector containing the error estimates}
  { iord   - the first k elements contain pointers to the error estimates}
  { nrmax  - maxerr = iord(nrmax)}
var
  errmax, errmin: extended;
  i, isucc, j, jbnd, jupbn, k : integer;
label
  done;
begin
  {Ref: Quadpack[34], subroutine dqpsrt in dqpsrt.f}

  if last<=2 then begin
    iord[1] := 1;
    iord[2] := 2;
    goto done;
  end;

  {This part of the routine is only executed if, due to a difficult integrand,}
  {subdivision increased the error estimate. in the normal case the insert    }
  {procedure should start after the nrmax-th largest error estimate.          }
  errmax := elist[maxerr];
  while (nrmax>1) and (errmax > elist[iord[nrmax-1]]) do begin
    iord[nrmax] := iord[nrmax-1];
    dec(nrmax);
  end;

  {Compute the number of elements in the list to be maintained in descending}
  {order. This number depends on the number of subdivisions still allowed.  }
  jupbn := last;
  if last > 2 + limit div 2 then jupbn := limit+3-last;
  errmin := elist[last];
  {Insert errmax by traversing the list top-down, starting}
  {comparison from the element elist[iord[nrmax+1]]       }
  jbnd := jupbn-1;
  for i:=nrmax+1 to jbnd do begin
    isucc := iord[i];
    if errmax > elist[isucc] then begin
     {insert errmin by traversing the list bottom-up}
     iord[i-1] := maxerr;
     k := jbnd;
     for j:=i to jbnd do begin
       isucc := iord[k];
       if errmin < elist[isucc] then begin
         iord[k+1] := last;
         goto done;
       end;
       iord[k+1] := isucc;
       dec(k);
     end;
     iord[i] := last;
     goto done;
    end;
    iord[i-1] := isucc;
  end;
  iord[jbnd]  := maxerr;
  iord[jupbn] := last;

done:
  {set maxerr and ermax}
  maxerr := iord[nrmax];
  ermax  := elist[maxerr];
end;


type
  TEpsArray = array[1..52] of extended;
  TEps3Res  = array[1..3]  of extended;


{---------------------------------------------------------------------------}
procedure qelg(var n: integer;        {index of the new element in the 1. column of the epsilon table}
               var epstab: TEpsArray; {the two lower diagonals of the triangular epsilon table}
               var result: extended;  {resulting approximation to the integral}
               var abserr: extended;  {estimate of the absolute error}
               var res3la: TEps3Res;  {vector containing the last 3 results}
               var nres: integer);    {number of calls to the routine (should be zero at first call)}
  {-Determine the limit of a given sequence of approximations by means of the epsilon algorithm}
const
  limexp = 50; {maximum number of elements the epsilon table can contain}
var
  e0, e1, e2, e3: extended; {the 4 elements on which the computation of a}
                            {new element in the epsilon table is based}
  error: extended;          {abs(e1-e0)+abs(e2-e1)+abs(new-e2)}
  delta1, delta2, delta3: extended;
  err1, err2, err3: extended;
  tol1, tol2, tol3: extended;
  e1abs, res, ss: extended;
  i, ib, ib2, k1, k2, k3, num: integer;
  newelm: integer; {number of elements to be computed in the new diagonal}
label
  l50;
begin
  {Ref: Quadpack[34], subroutine dqelg in dqelg.f}

  inc(nres);
  abserr := MaxExtended;
  result := epstab[n];
  if n<3 then exit;

  epstab[n+2] := epstab[n];
  epstab[n] := MaxExtended;
  num := n;
  k1  := n;
  newelm := (n-1) div 2;
  for i:=1 to newelm do begin
    k2  := k1-1;
    k3  := k1-2;
    res := epstab[k1+2];
    e0  := epstab[k3];
    e1  := epstab[k2];
    e2  := res;
    e1abs  := abs(e1);
    delta2 := e2-e1;
    err2   := abs(delta2);
    tol2   := maxx(e1abs,abs(e2))*eps_x;
    delta3 := e1-e0;
    err3   := abs(delta3);
    tol3   := maxx(e1abs, abs(e0))*eps_x;
    if (err2<=tol2) and (err3<=tol3) then begin
      {if e0,e1,e2 are equal to within machine accuracy, convergence is assumed}
      result := res;
      abserr := maxx(err2+err3, 5.0*eps_x*abs(result));
      exit;
    end;

    e3 := epstab[k1];
    epstab[k1] := e1;

    delta1 := e1-e3;
    err1   := abs(delta1);
    tol1   := maxx(e1abs, abs(e3))*eps_x;

    {if two elements are very close to each other, omit}
    {a part of the table by adjusting the value of n}
    if (err1<=tol1) or (err2<=tol2) or (err3<=tol3) then begin
      n := 2*i-1;
      goto l50;
    end;
    ss := 1.0/delta1 + 1.0/delta2 - 1.0/delta3;
    {Test to detect irregular behaviour in the table, and eventually}
    {omit a part of the table adjusting the value of n.}
    if abs(ss*e1) < 1e-4 then begin
      n := 2*i-1;
      goto l50;
    end;
    res := e1+1.0/ss;
    epstab[k1] := res;
    k1 := k1-2;
    error := err2 + abs(res-e2) + err3;
    if error <= abserr then begin
      abserr := error;
      result := res
    end;
  end;

l50:
  {shift the table}
  if n=limexp then n := 2*(limexp div 2)-1;
  if odd(num) then ib := 1 else ib := 2;

  for i:=1 to newelm+1 do begin
    ib2 := ib+2;
    epstab[ib] := epstab[ib2];
    ib := ib2
  end;
  for i:=1 to n do epstab[i] := epstab[num-n+i];
  if nres<4 then begin
    res3la[nres] := result;
    abserr := MaxExtended;
  end
  else begin
    {compute error estimate}
    abserr    := abs(result-res3la[3]) + abs(result-res3la[2]) + abs(result-res3la[1]);
    res3la[1] := res3la[2];
    res3la[2] := res3la[3];
    res3la[3] := result;
    abserr    := maxx(abserr, 5.0*eps_x*abs(result));
  end
end;


{---------------------------------------------------------------------------}
procedure qagse(f: TFuncX; a, b, epsabs, epsrel: extended;
                limit: integer; var result, abserr: extended;
                var neval: longint; var ier: integer;
                var alist, blist, elist, rlist: TQXLimArr;
                var iord: TQILimArr; var last: integer);
  {-Global adaptive quadrature of f over (a,b) based on 21-point Gauss-Kronrod}
  { rule for the subintervals, with acceleration by Wynn's epsilon algorithm.}

  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the number of subintervals              }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  { alist  - left end points of the subintervals                   }
  { blist  - right end points of the subintervals                  }
  { rlist  - integral approximations on the subintervals           }
  { elist  - absolute error estimates on the subintervals          }
  { iord   - first elements contain pointers to the error estimates}
  { last   - number of subintervals actually produced              }

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine the estimates for integral  }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges. If possible, }
  {        an appropriate special-purpose integrator should be used, which }
  {        is designed for handling the type of difficulty involved.       }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved. The error may be       }
  {        under-estimated.                                                }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 4 The algorithm does not converge. Roundoff error is detected in  }
  {        the extrapolation table. It is presumed that the requested      }
  {        tolerance cannot be achieved, and that the returned result is   }
  {        the best which can be obtained.                                 }
  {                                                                        }
  {    = 5 The integral is probably divergent, or slowly convergent. It    }
  {        must be noted that divergence can occur with any other value of }
  {        ier.                                                            }
  {                                                                        }
  {    = 6 The input is invalid, because epsabs <= 0 and epsrel < 50*eps_x.}
  {        result, abserr, neval, last are set to zero.                    }

var
  abseps,            {estimate of absolute error from epstab}
  area,              {sum of the integrals over the subintervals}
  erlarg,            {sum of the errors over the intervals larger than the smallest interval considered up to now}
  erlast,            {error on the interval currently subdivided}
  errbnd,            {requested accuracy max(epsabs,epsrel*abs(result))}
  errmax,            {elist[maxerr]}
  errsum,            {sum of the errors over the subintervals}
  small: extended;   {length of the smallest interval considered up to now, multiplied by 1.5}
  area1, area12, area2, a1, a2, b1, b2, correc, defabs,
  defab1, defab2, dres, ertest, resabs,reseps,
  error1, error2, erro12: extended;
  maxerr,            {pointer to the interval with largest error estimate}
  nres,              {number of calls to the extrapolation routine}
  nrmax,
  numrl2: integer;   {number of elements currently in rlist2}
  ktmin, id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn: integer;
  extrap,            {routine is attempting to perform extrapolation}
  noext: boolean;    {extrapolation is no longer allowed (true value)}
  res3la: TEps3Res;
  rlist2: TEpsArray; {the part of the epsilon table which is still needed for further computations}

label
  l90,l100,l115,l130;

begin
  {Ref: Quadpack[34], subroutine qagse in qagse.f}

  {test on validity of parameters}
  last   := 0;
  result := 0.0;
  abserr := 0.0;
  if (epsabs <= 0.0) and (epsrel < 50.0*eps_x) then begin
    ier := 6;
    neval := 0;
    exit;
  end;

  {first approximation to the integral}
  qk21(f, a, b, result, abserr, defabs, resabs);
  ier   := 0;
  ierro := 0;
  neval := 21;

  {test on accuracy}
  dres := abs(result);
  errbnd := maxx(epsabs, epsrel*dres);
  alist[1] := a;
  blist[1] := b;
  rlist[1] := result;
  elist[1] := abserr;
  iord[1] := 1;
  if (abserr <= 100.0*eps_x*defabs) and (abserr>errbnd) then ier := 2;
  if limit=1 then ier := 1;
  if (ier<>0) or ((abserr <= errbnd) and (abserr <> resabs)) or (abserr=0.0) then exit;

  {initialization}
  rlist2[1] := result;
  errmax := abserr;
  maxerr := 1;
  area   := result;
  errsum := abserr;
  erlarg := errsum;
  abserr := MaxExtended;
  small  := abs(b-a)*0.375;
  nrmax  := 1;
  nres   := 0;
  numrl2 := 2;
  ktmin  := 0;
  extrap := false;
  noext  := false;
  iroff1 := 0;
  iroff2 := 0;
  iroff3 := 0;
  ksgn := -1;
  if dres >= (1.0-50.0*eps_x)*defabs then ksgn := 1;

{$ifdef FPC}
  {Suppress warnings: Local variable does not seem to be initialized}
  correc := 0.0;
  ertest := 0.0;
{$endif}

  {main loop}
  last := 1;
  while last < limit do begin
    inc(last);
    {bisect the subinterval with the nrmax-th largest error estimate}
    a1 := alist[maxerr];
    b1 := 0.5*(alist[maxerr] + blist[maxerr]);
    a2 := b1;
    b2 := blist[maxerr];
    erlast := errmax;
    qk21(f, a1, b1, area1, error1, resabs, defab1);
    qk21(f, a2, b2, area2, error2, resabs, defab2);
    {improve previous approximations to integral and error and test for accuracy}
    area12 := area1  + area2;
    erro12 := error1 + error2;
    errsum := errsum + erro12 - errmax;
    area   := area + area12 - rlist[maxerr];
    if (defab1 <> error1) and (defab2 <> error2) then begin
      if (abs(rlist[maxerr]-area12) <= 1e-5*abs(area12)) and (erro12 >= 0.99*errmax) then begin
        if extrap then inc(iroff2)
        else inc(iroff1)
      end;
      if (last > 10) and (erro12 > errmax) then inc(iroff3)
    end;

    rlist[maxerr] := area1;
    rlist[last] := area2;
    errbnd := maxx(epsabs, epsrel*abs(area));
    {test for roundoff error and eventually set error flag}
    if (iroff1+iroff2 >= 10) or (iroff3>=20) then ier := 2;
    if iroff2>=5 then ierro := 3;
    {set error flag in the case that the number of subintervals equals limit}
    if last=limit then ier := 1;
    {set error flag in the case of bad integrand}
    {behaviour at a point of the integration range}
    if maxx(abs(a1),abs(b2)) <= (1.0+100.0*eps_x)*(abs(a2)+1000.0*MinExtended) then ier := 4;

    {append the newly-created intervals to the list}
    if error2 <= error1 then begin
      alist[last]   := a2;
      blist[maxerr] := b1;
      blist[last]   := b2;
      elist[maxerr] := error1;
      elist[last]   := error2;
    end
    else begin
      alist[maxerr] := a2;
      alist[last]   := a1;
      blist[last]   := b1;
      rlist[maxerr] := area2;
      rlist[last]   := area1;
      elist[maxerr] := error2;
      elist[last]   := error1;
    end;
    {call subroutine dqpsrt to maintain the descending ordering}
    {in the list of error estimates and select the subinterval }
    {with nrmax-th largest error estimate (to be bisected next)}
    qpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);

    if errsum <= errbnd then goto l115;
    if ier<>0 then goto l100;

    if last=2 then begin
      erlarg := errsum;
      ertest := errbnd;
      rlist2[2] := area
    end
    else if not noext then begin
      erlarg := erlarg-erlast;
      if abs(b1-a1) > small then erlarg := erlarg+erro12;
      if not extrap then begin
        {test whether the interval to be bisected next is the smallest interval}
        if abs(blist[maxerr]-alist[maxerr]) > small then goto l90;
        nrmax := 2;
        extrap := true;
      end;
      if (ierro<>3) and (erlarg>ertest) then begin
        {The smallest interval has the largest error.            }
        {Before bisecting decrease the sum of the errors over the}
        {larger intervals (erlarg) and perform extrapolation.    }
        id := nrmax;
        jupbnd := last;
        if last > 2 + limit div 2 then jupbnd := limit+3-last;
        for k:=id to jupbnd do begin
          maxerr := iord[nrmax];
          errmax := elist[maxerr];
          if abs(blist[maxerr]-alist[maxerr]) > small then goto l90;
          inc(nrmax);
        end;
      end;
      {perform extrapolation.}
      numrl2 := numrl2+1;
      rlist2[numrl2] := area;
      qelg(numrl2, rlist2, reseps, abseps, res3la, nres);
      ktmin := ktmin+1;
      if (ktmin > 5) and (abserr < 1e-3*errsum) then ier := 5;
      if abseps < abserr then begin
        ktmin  := 0;
        abserr := abseps;
        result := reseps;
        correc := erlarg;
        ertest := maxx(epsabs,epsrel*abs(reseps));
        if abserr <= ertest then goto l100;
      end;
      {prepare bisection of the smallest interval}
      if numrl2=1 then noext := true;
      if ier=5 then goto l100;
      maxerr := iord[1];
      errmax := elist[maxerr];
      nrmax  := 1;
      extrap := false;
      small  := 0.5*small;
      erlarg := errsum;
    end;
l90:
  end; {main loop}

l100:
  {set final result and error estimate.}
  if abserr=MaxExtended then goto l115;
  if ier+ierro<>0 then begin
    if ierro=3 then abserr := abserr+correc;
    if ier=0 then ier := 3;
    if (result=0.0) or (area=0.0) then begin
      if abserr>errsum then goto l115;
      if area=0.0      then goto l130;
    end;
  end;

  {test on divergence}
  if (ksgn=-1) and (maxx(abs(result),abs(area)) <= defabs*0.01) then goto l130;
  if (0.01 > result/area) or (result/area > 100.0) or (errsum>abs(area)) then ier := 6;
  goto l130;

l115:
  {compute global integral sum}
  result := 0.0;
  for k:=1 to last do result := result+rlist[k];
  abserr := errsum;

l130:
  if ier>2  then ier := ier-1;
  if last>0 then neval := 42*longint(last)-21;
end;


{---------------------------------------------------------------------------}
procedure qags(f: TFuncX; a, b, epsabs, epsrel: extended; limit: integer;
               var result, abserr: extended; var neval: longint; var ier: integer);
  {-Global adaptive quadrature of f over (a,b) based on 21-point Gauss-Kronrod}
  { rule for the subintervals, with acceleration by Wynn's epsilon algorithm.}
  { Simplified user interface to procedure qagse}

  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the no. of subintervals, 0: use DefLimit}
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine the estimates for integral  }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges. If possible, }
  {        an appropriate special-purpose integrator should be used, which }
  {        is designed for handling the type of difficulty involved.       }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved. The error may be       }
  {        under-estimated.                                                }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 4 The algorithm does not converge. Roundoff error is detected in  }
  {        the extrapolation table. It is presumed that the requested      }
  {        tolerance cannot be achieved, and that the returned result is   }
  {        the best which can be obtained.                                 }
  {                                                                        }
  {    = 5 The integral is probably divergent, or slowly convergent. It    }
  {        must be noted that divergence can occur with any other value of }
  {        ier.                                                            }
  {                                                                        }
  {    = 6 The input is invalid, because epsabs <= 0 and epsrel < 50*eps_x.}
  {        result, abserr, neval, last are set to zero.                    }
  {                                                                        }
  {    = 7 The input is invalid, limit < 0 or limit > QMAXLIM              }
  {        result, abserr, neval, last are set to zero.                    }
  {                                                                        }
  {    = 8 Dynamic list vectors cannot be allocated.                       }
  {        result, abserr, neval, last are set to zero.                    }
var
  last: integer;
  lsx : word;
  alist, blist, elist, rlist: ^TQXLimArr;
  iord: ^TQILimArr;
begin

  {set error return values}
  result := 0.0;
  abserr := 0.0;
  neval  := 0;

  if limit=0 then limit := DefLimit;
  if (limit<1) or (limit>QMAXLIM) then begin
    ier := 7;
    exit;
  end;
  lsx   := word(limit)*sizeof(extended);
  alist := malloc(lsx);
  blist := malloc(lsx);
  elist := malloc(lsx);
  rlist := malloc(lsx);
  iord  := malloc(limit*sizeof(integer));
  if (alist=nil) or (blist=nil) or (elist=nil) or (rlist=nil) or (iord=nil) then begin
    ier := 8;
  end
  else begin
    qagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval,ier,
          alist^, blist^, elist^, rlist^, iord^, last);
  end;
  mfree(alist, lsx);
  mfree(blist, lsx);
  mfree(elist, lsx);
  mfree(rlist, lsx);
  mfree(iord,  limit*sizeof(integer));
end;


{---------------------------------------------------------------------------}
procedure qagie(f: TFuncX; bound: extended; inf: integer;epsabs, epsrel: extended; limit: integer;
                var result, abserr: extended;
                var neval: longint;
                var ier: integer;
                var alist, blist, elist, rlist: TQXLimArr;
                var iord: TQILimArr; var last: integer);
  {-Global adaptive quadrature of f over an infinite interval based on}
  { transformed 15-point Gauss-Kronrod rule for the subintervals, with}
  { acceleration by Wynn's epsilon algorithm.}

  { f      - function defining the integrand                       }
  { bound  - finite bound of integration ran                       }
  { inf    - indicating the kind of integration range involved     }
  {              inf =  1 corresponds to  (bound, +infinity)       }
  {              inf = -1             to  (-infinity, bound)       }
  {              inf =  2             to  (-infinity, +infinity)   }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the number of subintervals, must be > 1 }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  { alist  - left end points of the subintervals                   }
  { blist  - right end points of the subintervals                  }
  { rlist  - integral approximations on the subintervals           }
  { elist  - absolute error estimates on the subintervals          }
  { iord   - first elements contain pointers to the error estimates}
  { last   - number of subintervals actually produced              }

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine the estimates for integral  }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges. If possible, }
  {        an appropriate special-purpose integrator should be used, which }
  {        is designed for handling the type of difficulty involved.       }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved. The error may be       }
  {        under-estimated.                                                }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 4 The algorithm does not converge. Roundoff error is detected in  }
  {        the extrapolation table. It is presumed that the requested      }
  {        tolerance cannot be achieved, and that the returned result is   }
  {        the best which can be obtained.                                 }
  {                                                                        }
  {    = 5 The integral is probably divergent, or slowly convergent. It    }
  {        must be noted that divergence can occur with any other value of }
  {        ier.                                                            }
  {                                                                        }
  {    = 6 The input is invalid, because epsabs <= 0 and epsrel < 50*eps_x.}
  {        result, abserr, neval, last are set to zero.                    }

var
  abseps,            {estimate of absolute error from epstab}
  area,              {sum of the integrals over the subintervals}
  erlarg,            {sum of the errors over the intervals larger than the smallest interval considered up to now}
  erlast,            {error on the interval currently subdivided}
  errbnd,            {requested accuracy max(epsabs,epsrel*abs(result))}
  errmax,            {elist[maxerr]}
  errsum,            {sum of the errors over the subintervals}
  small: extended;   {length of the smallest interval considered up to now, multiplied by 1.5}
  area1, area12, area2, a1, a2, b1, b2, correc, defabs,
  defab1, defab2, dres, ertest, resabs,reseps,
  error1, error2, erro12: extended;
  maxerr,            {pointer to the interval with largest error estimate}
  nres,              {number of calls to the extrapolation routine}
  nrmax,
  numrl2: integer;   {number of elements currently in rlist2}
  ktmin, id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn: integer;
  extrap,            {routine is attempting to perform extrapolation}
  noext: boolean;    {extrapolation is no longer allowed (true value)}
  res3la: TEps3Res;
  rlist2: TEpsArray; {the part of the epsilon table which is still needed for further computations}
label
  l90,l100,l115,l130;
begin
  {Ref: Quadpack[34], subroutine dqagie in dqagie.f}

  {test on validity of parameter}
  if (epsabs<=0.0) and (epsrel < 50.0*eps_x) then begin
    ier    := 6;
    neval  := 0;
    last   := 0;
    result := 0.0;
    abserr := 0.0;
    exit;
  end;

  {first approximation to the integral}
  {-----------------------------------}

  {determine the interval to be mapped onto (0,1).}
  {if inf=2 the integral is computed as i = i1+i2, where}
  {  i1 = integral of f over (-infinity,0), }
  {  i2 = integral of f over (0,+infinity). }

  if (inf=2) then bound := 0;
  qk15i(f, bound, inf, 0, 1, result, abserr, defabs, resabs);

  {test on accuracy}
  dres := abs(result);
  errbnd := maxx(epsabs,epsrel*dres);

  ier := 0;
  if (abserr<=100.0*eps_x*defabs) and (abserr>errbnd) then ier := 2;
  if limit=1 then ier := 1;
  if (ier<>0) or ((abserr<=errbnd) and (abserr<>resabs)) or (abserr=0) then exit;

  {initialization}
  alist[1] := 0;
  blist[1] := 1;
  rlist[1] := result;
  elist[1] := abserr;
  iord[1]  := 1;
  rlist2[1]:= result;
  errmax   := abserr;
  maxerr   := 1;
  area     := result;
  errsum   := abserr;
  erlarg   := errsum;
  abserr   := MaxExtended;
  small    := 0.375;
  nrmax    := 1;
  nres     := 0;
  ktmin    := 0;
  numrl2   := 2;
  extrap   := false;
  noext    := false;
  ierro    := 0;
  iroff1   := 0;
  iroff2   := 0;
  iroff3   := 0;
  if dres>=(1.0-50.0*eps_x)*defabs then ksgn := 1 else ksgn := -1;

{$ifdef FPC}
  {Suppress warnings: Local variable does not seem to be initialized}
  correc := 0.0;
  ertest := 0.0;
{$endif}

  last := 1;
  {main loop}
  while last<limit do begin
    {bisect the subinterval with nrmax-th largest error estimate}
    inc(last);
    a1 := alist[maxerr];
    b1 := 0.5*(alist[maxerr] + blist[maxerr]);
    a2 := b1;
    b2 := blist[maxerr];
    erlast := errmax;
    qk15i(f, bound, inf, a1, b1, area1, error1, resabs, defab1);
    qk15i(f, bound, inf, a2, b2, area2, error2, resabs, defab2);
    {improve previous approximations to integral and error and test for accuracy}
    area12 := area1+area2;
    erro12 := error1+error2;
    errsum := errsum+erro12-errmax;
    area   := area+area12-rlist[maxerr];
    if (defab1<>error1) and (defab2<>error2) then begin
      if (abs(rlist[maxerr]-area12)<=1e-5*abs(area12)) and (erro12>=0.99*errmax) then begin
        if extrap then inc(iroff2)
        else inc(iroff1);
      end;
      if (last>10) and (erro12>errmax) then inc(iroff3)
    end;
    rlist[maxerr] := area1;
    rlist[last] := area2;
    errbnd := maxx(epsabs,epsrel*abs(area));

    {test for roundoff error and eventually set error flag}
    if (iroff1+iroff2>=10) or (iroff3>=20) then ier := 2;
    if iroff2>=5  then ierro := 3;
    if last=limit then ier := 1;

    {set error flag in the case of bad integrand behaviour}
    {at some points of the integration range.}
    if (maxx(abs(a1),abs(b2)) <= (1.0+100.0*eps_x)*(abs(a2)+1000.0*MinExtended)) then ier := 4;

    {append the newly-created intervals to the list}
    if error2 <= error1 then begin
      alist[last]   := a2;
      blist[maxerr] := b1;
      blist[last]   := b2;
      elist[maxerr] := error1;
      elist[last]   := error2;
    end
    else begin
      alist[maxerr] := a2;
      alist[last]   := a1;
      blist[last]   := b1;
      rlist[maxerr] := area2;
      rlist[last]   := area1;
      elist[maxerr] := error2;
      elist[last]   := error1;
    end;
    {call subroutine dqpsrt to maintain the descending ordering}
    {in the list of error estimates and select the subinterval }
    {with nrmax-th largest error estimate (to be bisected next)}
    qpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);

    if errsum<=errbnd then goto l115;
    if ier<>0 then goto l100;

    if last=2 then begin
      erlarg := errsum;
      ertest := errbnd;
      rlist2[2] := area
    end
    else if not noext then begin
      erlarg := erlarg-erlast;
      if abs(b1-a1) > small then erlarg := erlarg+erro12;
      if not extrap then begin
        {test whether the interval to be bisected next is the smallest interval}
        if abs(blist[maxerr]-alist[maxerr]) > small then goto l90;
        nrmax := 2;
        extrap := true;
      end;
      if (ierro<>3) and (erlarg>ertest) then begin
        {The smallest interval has the largest error.            }
        {Before bisecting decrease the sum of the errors over the}
        {larger intervals (erlarg) and perform extrapolation.    }
        id := nrmax;
        jupbnd := last;
        if last > 2 + limit div 2 then jupbnd := limit+3-last;
        for k:=id to jupbnd do begin
          maxerr := iord[nrmax];
          errmax := elist[maxerr];
          if abs(blist[maxerr]-alist[maxerr]) > small then goto l90;
          inc(nrmax);
        end;
      end;
      {perform extrapolation.}
      numrl2 := numrl2+1;
      rlist2[numrl2] := area;
      qelg(numrl2, rlist2, reseps, abseps, res3la, nres);
      ktmin := ktmin+1;
      if (ktmin > 5) and (abserr < 1e-3*errsum) then ier := 5;
      if abseps < abserr then begin
        ktmin  := 0;
        abserr := abseps;
        result := reseps;
        correc := erlarg;
        ertest := maxx(epsabs,epsrel*abs(reseps));
        if abserr <= ertest then goto l100;
      end;
      {prepare bisection of the smallest interval}
      if numrl2=1 then noext := true;
      if ier=5 then goto l100;
      maxerr := iord[1];
      errmax := elist[maxerr];
      nrmax  := 1;
      extrap := false;
      small  := small/2;
      erlarg := errsum;
    end;
l90:
  end; {main loop}

l100:
  {set final result and error estimate}
  if abserr=MaxExtended then goto l115;
  if ier+ierro<>0 then begin
    if ierro=3 then abserr := abserr+correc;
    if ier=0 then ier := 3;
    if (result=0.0) or (area=0.0) then begin
      if abserr>errsum then goto l115;
      if area=0.0      then goto l130;
    end;
  end;

  {test on divergence}
  if (ksgn=-1) and (maxx(abs(result),abs(area)) <= defabs*0.01) then goto l130;
  if (0.01 > result/area) or (result/area > 100.0) or (errsum>abs(area)) then ier := 6;
  goto l130;

l115:
  {compute global integral sum}
  result := 0.0;
  for k:=1 to last do result := result+rlist[k];
  abserr := errsum;

l130:
  if ier>2  then ier := ier-1;
  if last>0 then neval := 30*longint(last)-15;
  if inf=2  then neval := 2*neval;
end;


{---------------------------------------------------------------------------}
procedure qagi(f: TFuncX; bound: extended;inf: integer; epsabs, epsrel: extended; limit: integer;
               var result, abserr: extended; var neval: longint; var ier: integer);
  {-Global adaptive quadrature of f over an infinite interval based on trans- }
  { formed 15-point Gauss-Kronrod rule for the subintervals, with acceleration}
  { by Wynn's epsilon algorithm. Simplified user interface to procedure qagie.}

  { f      - function defining the integrand                       }
  { bound  - finite bound of integration ran                       }
  { inf    - indicating the kind of integration range involved     }
  {              inf =  1 corresponds to  (bound, +infinity)       }
  {              inf = -1             to  (-infinity, bound)       }
  {              inf =  2             to  (-infinity, +infinity)   }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the no. of subintervals, 0: use DefLimit}
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine the estimates for integral  }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges. If possible, }
  {        an appropriate special-purpose integrator should be used, which }
  {        is designed for handling the type of difficulty involved.       }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved. The error may be       }
  {        under-estimated.                                                }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 4 The algorithm does not converge. Roundoff error is detected in  }
  {        the extrapolation table. It is presumed that the requested      }
  {        tolerance cannot be achieved, and that the returned result is   }
  {        the best which can be obtained.                                 }
  {                                                                        }
  {    = 5 The integral is probably divergent, or slowly convergent. It    }
  {        must be noted that divergence can occur with any other value of }
  {        ier.                                                            }
  {                                                                        }
  {    = 6 The input is invalid, because epsabs <= 0 and epsrel < 50*eps_x.}
  {        result, abserr, neval, last are set to zero.                    }
var
  last: integer;
  lsx : word;
  alist, blist, elist, rlist: ^TQXLimArr;
  iord: ^TQILimArr;
begin
  if limit=0 then limit := DefLimit;
  if (limit<2) or (limit>QMAXLIM) then begin
    ier := 7;
    exit;
  end;
  lsx   := word(limit)*sizeof(extended);
  alist := malloc(lsx);
  blist := malloc(lsx);
  elist := malloc(lsx);
  rlist := malloc(lsx);
  iord  := malloc(limit*sizeof(integer));
  if (alist=nil) or (blist=nil) or (elist=nil) or (rlist=nil) or (iord=nil) then begin
    ier := 8;
  end
  else begin
    qagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,neval,ier,
          alist^, blist^, elist^, rlist^, iord^, last);
  end;
  mfree(alist, lsx);
  mfree(blist, lsx);
  mfree(elist, lsx);
  mfree(rlist, lsx);
  mfree(iord,  limit*sizeof(integer));
end;


{---------------------------------------------------------------------------}
procedure quagk(f: TFuncX; a,b,epsabs: extended; var result,abserr: extended; var ier: integer);
  {-Global adaptive quadrature of f over (a,b) based on Gauss-Kronrod rules}
  { for the subintervals, with acceleration by Wynn's epsilon algorithm.}
  { Simplified user interface to procedure qags and qagi.}

  { f      - function defining the integrand                       }
  { a      - lower limit of integration, may be infinite           }
  { b      - upper limit of integration, may be infinite           }
  { epsabs - absolute accuracy requested                           }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { ier    - return code (success=0, error>0 see below)            }

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine the estimates for integral  }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges. If possible, }
  {        an appropriate special-purpose integrator should be used, which }
  {        is designed for handling the type of difficulty involved.       }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved. The error may be       }
  {        under-estimated.                                                }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 4 The algorithm does not converge. Roundoff error is detected in  }
  {        the extrapolation table. It is presumed that the requested      }
  {        tolerance cannot be achieved, and that the returned result is   }
  {        the best which can be obtained.                                 }
  {                                                                        }
  {    = 5 The integral is probably divergent, or slowly convergent. It    }
  {        must be noted that divergence can occur with any other value of }
  {        ier.                                                            }
  {                                                                        }
  {    = 6 The input is invalid, because epsabs <= 0 and epsrel < 50*eps_x.}
  {        result, abserr, last are set to zero.                           }
  {                                                                        }
  {    = 7 The input is invalid, limit < 0 or limit > QMAXLIM.             }
  {        result, abserr, last are set to zero.                           }
  {                                                                        }
  {    = 8 Dynamic list vectors cannot be allocated.                       }
  {        result, abserr, last are set to zero.                           }
  {                                                                        }
  {    = 9 At least one limit a or b is NaN, or infinite a=b.              }
  {        result, abserr, last are set to NaN_x.                          }
var
  epsrel: extended;
  neval: longint;
  infa,infb: boolean;
begin
  if epsabs<0.0 then epsabs := 0.0;
  if epsabs=0.0 then epsrel := 50.0*eps_x else epsrel := 0.0;
  infa := IsInf(a);
  infb := IsInf(b);
  if IsNan(a) or IsNan(b) or (infa and (a=b)) then begin
    ier := 9;
    result := NaN_x;
    abserr := NaN_x;
    exit;
  end;
  if infa and infb then begin
    {a and b are infinite but a <> b}
    qagi(f, 0, 2, epsabs, epsrel, 0, result, abserr, neval, ier);
    if a>b then result := -result;
  end
  else if infa then begin
    {a is infinite, b is finite}
    qagi(f, b, isign(a), epsabs, epsrel, 0, result, abserr, neval, ier);
    if a>0.0 then result := -result;
  end
  else if infb then begin
    {a is finite, b is infinite}
    qagi(f, a, isign(b), epsabs, epsrel, 0, result, abserr, neval, ier);
    if b<0.0 then result := -result;
  end
  else begin
    {a and b are finite}
    qags(f, a, b, epsabs, epsrel, 0, result, abserr, neval, ier);
  end;
end;


{---------------------------------------------------------------------------}
procedure qk15c(f: TFuncX; a,b,c: extended; var result, abserr, resabs, resasc: extended);
  {-Integration of f(x)/(x-c) over (a,b) with 15-point Gauss-Kronrod rule}
const
  wg : array[1..4] of extended = ( {weights of the 7-point Gauss rule}
         0.129484966168869693270611432679082,
         0.279705391489276667901467771423780,
         0.381830050505118944950369775488975,
         0.417959183673469387755102040816327);
  xgk: array[1..8] of extended = ( {abscissae of the 15-point Kronrod rule}
         0.991455371120812639206854697526329,
         0.949107912342758524526189684047851,
         0.864864423359769072789712788640926,
         0.741531185599394439863864773280788,
         0.586087235467691130294144838258730,
         0.405845151377397166906606412076961,
         0.207784955007898467600689403773245,
         0.0);
  wgk: array[1..8] of extended = ( {weights of the 15-point Kronrod rule}
         0.022935322010529224963732008058970,
         0.063092092629978553290700663189204,
         0.104790010322250183839876322541518,
         0.140653259715525918745189590510238,
         0.169004726639267902826583426598550,
         0.190350578064785409913256402421014,
         0.204432940075298892414161999234649,
         0.209482141084727828012999174891714);
var
  centr: extended;  {mid point of the interval}
  absc,
  absc1,
  absc2: extended;  {abscissa}
  fval1,
  fval2: extended;  {function value}
  hlgth: extended;  {half-length of the interval}
  resg:  extended;  {result of the 7-point Gauss formula}
  resk:  extended;  {result of the 15-point Kronrod formula}
  reskh: extended;  {approximation to the mean value of the transformed integrand}
  fc, fsum, t: extended;
  fv1, fv2: array[1..7] of extended;
  j,k: integer;
begin

  {Ref: Quadpack[34], subroutines dqk15 and dqk15w. Quadpack uses dqk15w and}
  {a special parameter passing mechanism with another separate function code}
  {derived from f. Here I simply use the normal routine dqk15 for f(x)/(x-c)}

  centr := 0.5*(a+b);
  hlgth := 0.5*(b-a);

  {compute the 15-point Kronrod approximation to the integral, and estimate the error}
  fc   := f(centr)/(centr-c);
  resg := wg[4]*fc;
  resk := wgk[8]*fc;
  resabs := abs(resk);
  for j:=1 to 3 do begin
    k := j*2;
    absc   := hlgth*xgk[k];
    absc1  := centr - absc;
    absc2  := centr + absc;
    fval1  := f(absc1)/(absc1-c);
    fval2  := f(absc2)/(absc2-c);
    fv1[k] := fval1;
    fv2[k] := fval2;
    fsum   := fval1  + fval2;
    resg   := resg   + wg[j]*fsum;
    resk   := resk   + wgk[k]*fsum;
    resabs := resabs + wgk[k]*(abs(fval1) + abs(fval2));
  end;
  for j:=1 to 4 do begin
    k := j*2-1;
    absc   := hlgth*xgk[k];
    absc1  := centr - absc;
    absc2  := centr + absc;
    fval1  := f(absc1)/(absc1-c);
    fval2  := f(absc2)/(absc2-c);
    fv1[k] := fval1;
    fv2[k] := fval2;
    fsum   := fval1  + fval2;
    resk   := resk   + wgk[k]*fsum;
    resabs := resabs + wgk[k]*(abs(fval1) + abs(fval2));
  end;

  reskh := 0.5*resk;
  resasc := wgk[8]*abs(fc-reskh);
  for j:=1 to 7 do begin
    resasc := resasc + wgk[j]*(abs(fv1[j]-reskh) + abs(fv2[j]-reskh));
  end;

  t := abs(hlgth);
  result := resk*hlgth;
  resabs := resabs*t;
  resasc := resasc*t;
  abserr := abs((resk-resg)*hlgth);

  {scale error}
  if (resasc<>0.0) and (abserr<>0.0) then begin
    t := sqrt(abs(abserr/resasc));
    abserr := resasc*minx(1.0,200.0*t*t*t)
  end;
  t := 50.0*eps_x;
  if resabs > MinExtended/t then abserr := maxx(t*resabs,abserr);
end;


{---------------------------------------------------------------------------}
procedure qc25gcc(f: TFuncX; a,b,c: extended; var result,abserr: extended);
  {-Integration of f(x)/(x-c) over (a,b) using the generalized Clenshaw-Curtis method}
const
  x: array[1..11] of extended = ( {x contains the values cos(k*pi/24)}
       0.991444861373810411144557526928563,
       0.965925826289068286749743199728897,
       0.923879532511286756128183189396788,
       0.866025403784438646763723170752936,
       0.793353340291235164579776961501299,
       0.707106781186547524400844362104849,
       0.608761429008720639416097542898164,
       0.500000000000000000000000000000000,
       0.382683432365089771728459984030399,
       0.258819045102520762348898837624048,
       0.130526192220051591548406227895489);
var
  v,cheb12: array[1..13] of extended;
  fval,cheb24: array[1..25] of extended;
var
  ak22,amom0,amom1,amom2,cc,centr,hlgth: extended;
  res12,res24,u: extended;
  alam,alam1,alam2,part1,part2,part3: extended;
var
  i,j,k: integer;
begin

  {25-point Clenshaw-Curtis rule for computating a Cauchy principal value }

  {Ref: Quadpack[34], subroutines dqc25c and dqcheb. The original Quadpack}
  {dqc25c uses dqcheb and dqk15w as subroutines. My routine has the code  }
  {from dqcheb inline and the dqk15w equivalent is called from the driver }
  {procedure, this minimizes stack space and is easier to read and debug. }

  cc := (2.0*c-b-a)/(b-a);
  hlgth := 0.5*(b-a);
  centr := 0.5*(b+a);
  fval[1] := 0.5*f(hlgth+centr);
  fval[13] := f(centr);
  fval[25] := 0.5*f(centr-hlgth);
  for i:=2 to 12 do begin
    u := hlgth*x[i-1];
    j := 26-i;
    fval[i] := f(u+centr);
    fval[j] := f(centr-u);
  end;

  {----------------------------------------------------------------}
  {Compute the Chebyshev series expansion, this is the inlined code}
  {from dqcheb.f for the function call dqcheb(x,fval,cheb12,cheb24)}
  for i:=1 to 12 do begin
    j:=26-i;
    v[i]:=fval[i]-fval[j];
    fval[i]:=fval[i]+fval[j];
  end;

  alam1 := v[1]-v[9];
  alam2 := x[6]*(v[3]-v[7]-v[11]);
  cheb12[4] := alam1+alam2;
  cheb12[10]:= alam1-alam2;

  alam1 := v[2] - v[8] - v[10];
  alam2 := v[4] - v[6] - v[12];
  alam  := x[3]*alam1 + x[9]*alam2;
  cheb24[4] := cheb12[4] + alam;
  cheb24[22]:= cheb12[4] - alam;

  alam := x[9]*alam1-x[3]*alam2;
  cheb24[10] := cheb12[10] + alam;
  cheb24[16] := cheb12[10] - alam;

  part1 := x[4]*v[5];
  part2 := x[8]*v[9];
  part3 := x[6]*v[7];
  alam1 := v[1] + part1 + part2;
  alam2 := x[2]*v[3] + part3 + x[10]*v[11];
  cheb12[2] := alam1 + alam2;
  cheb12[12]:= alam1 - alam2;

  alam := x[1]*v[2] + x[3]*v[4] + x[5]*v[6] + x[7]*v[8] + x[9]*v[10] + x[11]*v[12];
  cheb24[2]  := cheb12[2] + alam;
  cheb24[24] := cheb12[2] - alam;

  alam := x[11]*v[2] - x[9]*v[4] + x[7]*v[6] - x[5]*v[8] + x[3]*v[10] - x[1]*v[12];
  cheb24[12] := cheb12[12] + alam;
  cheb24[14] := cheb12[12] - alam;

  alam1 := v[1] - part1 + part2;
  alam2 := x[10]*v[3] - part3 + x[2]*v[11];
  cheb12[6] := alam1 + alam2;
  cheb12[8] := alam1 - alam2;

  alam := x[5]*v[2] - x[9]*v[4] - x[1]*v[6] - x[11]*v[8] + x[3]*v[10] + x[7]*v[12];
  cheb24[6]  := cheb12[6] + alam;
  cheb24[20] := cheb12[6] - alam;

  alam := x[7]*v[2] - x[3]*v[4] - x[11]*v[6] + x[1]*v[8] - x[9]*v[10] - x[5]*v[12];
  cheb24[8]  := cheb12[8] + alam;
  cheb24[18] := cheb12[8] - alam;

  for i:=1 to 6 do begin
    j := 14-i;
    v[i]    := fval[i] - fval[j];
    fval[i] := fval[i] + fval[j];
  end;

  alam1 := v[1] + x[8]*v[5];
  alam2 := x[4]*v[3];
  cheb12[3]  := alam1 + alam2;
  cheb12[11] := alam1 - alam2;
  cheb12[7]  := v[1]  - v[5];

  alam := x[2]*v[2] + x[6]*v[4] + x[10]*v[6];
  cheb24[3]  := cheb12[3] + alam;
  cheb24[23] := cheb12[3] - alam;

  alam := x[6]*(v[2] - v[4] - v[6]);
  cheb24[7]  := cheb12[7] + alam;
  cheb24[19] := cheb12[7] - alam;

  alam := x[10]*v[2] - x[6]*v[4] + x[2]*v[6];
  cheb24[11] := cheb12[11] + alam;
  cheb24[15] := cheb12[11] - alam;
  for i:=1 to 3 do begin
    j := 8-i;
    v[i]    := fval[i] - fval[j];
    fval[i] := fval[i] + fval[j];
  end;

  cheb12[5] := v[1] + x[8]*v[3];
  cheb12[9] := fval[1] - x[8]*fval[3];

  alam := x[4]*v[2];
  cheb24[5]  := cheb12[5] + alam;
  cheb24[21] := cheb12[5] - alam;

  alam := x[8]*fval[2] - fval[4];
  cheb24[9]  := cheb12[9] + alam;
  cheb24[17] := cheb12[9] - alam;
  cheb12[1]  := fval[1] + fval[3];

  alam := fval[2] + fval[4];
  cheb24[1]  := cheb12[1] + alam;
  cheb24[25] := cheb12[1] - alam;
  cheb12[13] := v[1] - v[3];
  cheb24[13] := cheb12[13];

  for i:=2 to 12 do cheb12[i] := cheb12[i]/6.0;
  cheb12[1]:=cheb12[1]/12.0;
  cheb12[13]:=cheb12[13]/12.0;

  for i:=2 to 24 do cheb24[i] := cheb24[i]/12.0;
  cheb24[1]  := cheb24[1]/24.0;
  cheb24[25] := cheb24[25]/24.0;
  {------------------- end of inlined code ------------------------}

  {the modified chebyshev moments are computed by forward}
  {recursion, using amom0 and amom1 as starting values.}
  amom0 := ln(abs((1.0-cc)/(1.0+cc)));
  amom1 := 2.0 + cc*amom0;
  res12 := cheb12[1]*amom0 + cheb12[2]*amom1;
  res24 := cheb24[1]*amom0 + cheb24[2]*amom1;
  for k:=3 to 13 do begin
    amom2 := 2.0*cc*amom1 - amom0;
    ak22  := sqr(k-2);
    if k and 1 = 0 then amom2 := amom2 - 4.0/(ak22-1.0);
    res12 := res12 + cheb12[k]*amom2;
    res24 := res24 + cheb24[k]*amom2;
    amom0 := amom1;
    amom1 := amom2;
  end;
  for k:=14 to 25 do begin
    amom2 := 2.0*cc*amom1 - amom0;
    ak22  := sqr(k-2);
    if k and 1 = 0 then amom2 := amom2 - 4.0/(ak22-1.0);
    res24 := res24 + cheb24[k]*amom2;
    amom0 := amom1;
    amom1 := amom2;
  end;
  result := res24;
  abserr := abs(res24-res12);
end;


{---------------------------------------------------------------------------}
procedure qc25c(f: TFuncX; a,b,c: extended; var result,abserr: extended; var krul,neval: integer);
  {-Integration of f(x)/(x-c) over (a,b). Selects between 15-point}
  { Gauss-Kronrod rule and the generalized Clenshaw-Curtis method.}
  { krul is a flag which is decreased by 1 if the 15-point Gauss-Kronrod}
  { scheme has been used and is reliable}
var
  resabs,resasc: extended;
begin
  if abs(2.0*c-b-a) >= 1.1*abs(b-a) then begin
    {apply the 15-point gauss-kronrod scheme}
    qk15c(f,a,b,c,result,abserr,resabs,resasc);
    neval := 15;
    if resasc<>abserr then krul := krul-1;
  end
  else begin
    {use the generalized clenshaw-curtis method}
    qc25gcc(f,a,b,c,result,abserr);
    neval := 25;
  end;
end;


{---------------------------------------------------------------------------}
procedure qawce(f: TFuncX; a, b, c, epsabs, epsrel: extended;
                limit: integer; var result, abserr: extended;
                var neval: longint; var ier: integer;
                var alist, blist, elist, rlist: TQXLimArr;
                var iord: TQILimArr; var last: integer);
  {-Adaptive quadrature of the function f(x)/(x-c) over the finite interval}
  { (a,b) with the singularity at c with c<>a, c<>b. The routine calculates}
  { an approximation result to the Cauchy principal value. Parameters:}

  {#F}
  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { c      - singularity, ier=6 if c=a or c=b                      }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the number of subintervals              }
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  { alist  - left end points of the subintervals                   }
  { blist  - right end points of the subintervals                  }
  { rlist  - integral approximations on the subintervals           }
  { elist  - absolute error estimates on the subintervals          }
  { iord   - first elements contain pointers to the error estimates}
  { last   - number of subintervals actually produced              }
  {#F}

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine, the estimates for integral }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges.              }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved.                        }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 6 The input is invalid, because c = a or c = b, or epsabs <= 0 and}
  {        epsrel < 50*eps_x. result, abserr, neval, last are set to zero. }

var
  aa,bb,a1,a2,b1,b2,area,area1,area2,area12,
  error1,error2,erro12,errbnd,errmax,errsum: extended;
  maxerr,nrmax,iroff1,iroff2,k: integer;
  krule,nev: integer;
begin
  {test on validity of parameters}
  last := 0;
  result := 0.0;
  abserr := 0.0;
  if (c=a) or (c=b) or ((epsabs<=0.0) and (epsrel < 50.0*eps_x)) then begin
    ier := 6;
    neval := 0;
    exit;
  end;

  {first approximation to the integral}
  if a<=b then begin
    aa := a;
    bb := b;
  end
  else begin
    aa := b;
    bb := a;
  end;
  ier := 0;

  krule := 1;
  qc25c(f,aa,bb,c,result,abserr,krule,nev);
  neval := nev;
  last  := 1;
  rlist[1] := result;
  elist[1] := abserr;
  iord[1]  := 1;
  alist[1] := a;
  blist[1] := b;

  {test on accuracy}
  errbnd := maxx(epsabs,epsrel*abs(result));
  if limit=1 then ier := 1;
  if (abserr<minx(0.01*abs(result),errbnd)) or (ier=1) then begin
    if aa=b then result := -result;
    exit;
  end;

  {initialization}
  alist[1] := aa;
  blist[1] := bb;
  rlist[1] := result;
  errmax := abserr;
  maxerr := 1;
  area := result;
  errsum := abserr;
  nrmax := 1;
  iroff1 := 0;
  iroff2 := 0;

  {main loop}
  last := 1;
  repeat
    inc(last);

    {bisect the subinterval with nrmax-th largest error estimate}
    a1 := alist[maxerr];
    b1 := 0.5*(alist[maxerr] + blist[maxerr]);
    b2 := blist[maxerr];
    if (c<=b1) and (c>a1) then b1 := 0.5*(c+b2);
    if (c >b1) and (c<b2) then b1 := 0.5*(a1+c);
    a2 := b1;
    krule := 2;
    qc25c(f,a1,b1,c,area1,error1,krule,nev);
    neval := neval+nev;
    qc25c(f,a2,b2,c,area2,error2,krule,nev);
    neval := neval+nev;

    {improve previous approximations to integral and error and test for accuracy}
    area12 := area1  + area2;
    erro12 := error1 + error2;
    errsum := errsum + erro12 - errmax;
    area   := area   + area12 - rlist[maxerr];

    if (abs(rlist[maxerr]-area12) < 1e-5*abs(area12)) and (erro12 >= 0.99*errmax) and (krule=0) then inc(iroff1);
    if (last>10) and (erro12 > errmax) and (krule=0) then inc(iroff2);
    rlist[maxerr] := area1;
    rlist[last]   := area2;
    errbnd := maxx(epsabs,epsrel*abs(area));

    if errsum > errbnd then begin
      {test for roundoff error and eventually set error flag}
      if (iroff1 >= 6) and (iroff2 > 20) then ier := 2;
      {set error flag in the case that number of interval bisections exceeds limit}
      if last=limit then ier := 1;
      {set error flag in the case of bad integrand}
      {behaviour at a point of the integration range}
      if maxx(abs(a1),abs(b2)) <= (1.0+100.0*eps_x)*(abs(a2)+1000.0*MinExtended) then ier := 3;
    end;
    {append the newly-created intervals to the list}
    if error2 <= error1 then begin
      alist[last]   := a2;
      blist[maxerr] := b1;
      blist[last]   := b2;
      elist[maxerr] := error1;
      elist[last]   := error2;
    end
    else begin
      alist[maxerr] := a2;
      alist[last]   := a1;
      blist[last]   := b1;
      rlist[maxerr] := area2;
      rlist[last]   := area1;
      elist[maxerr] := error2;
      elist[last]   := error1;
    end;
    {call subroutine  qpsrt to maintain the descending ordering}
    {in the list of error estimates and select the subinterval }
    {with nrmax-th largest error estimate (to be bisected next)}
    qpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);

  until (ier<>0) or (errsum <= errbnd) or (last>=limit);

  {compute final result}
  result := 0.0;
  for k:=1 to last do result := result + rlist[k];
  abserr := errsum;
  if aa=b then result := -result;
end;


{---------------------------------------------------------------------------}
procedure qawc(f: TFuncX; a, b, c, epsabs, epsrel: extended;
               limit: integer; var result, abserr: extended;
               var neval: longint; var ier: integer);
  {-Adaptive quadrature of the function f(x)/(x-c) over the finite interval}
  { (a,b) with the singularity at c with c<>a, c<>b. The routine calculates}
  { an approximation result to the Cauchy principal value.  Simplified user}
  { interface to procedure qawce. Parameters:}

  {#F}
  { f      - function defining the integrand                       }
  { a      - lower limit of integration                            }
  { b      - upper limit of integration                            }
  { c      - singularity, ier=6 if c=a or c=b                      }
  { epsabs - absolute accuracy requested                           }
  { epsrel - relative accuracy requested                           }
  { limit  - upperbound on the no. of subintervals, 0: use DefLimit}
  { result - approximation to the integral                         }
  { abserr - estimate of the modulus of the absolute error         }
  { neval  - number of integrand evaluations                       }
  { ier    - return code (success=0, error>0 see below)            }
  {#F}

  {ier = 0 Normal and reliable termination of the routine. It is assumed   }
  {        that the requested accuracy has been achieved.                  }
  {                                                                        }
  {ier > 0 Abnormal termination of the routine, the estimates for integral }
  {        and error are less reliable. It is assumed that the requested   }
  {        accuracy has not been achieved.                                 }
  {                                                                        }
  {    = 1 Maximum number of subdivisions allowed has been achieved. One   }
  {        can allow more subdivisions by increasing the value of limit    }
  {        (and taking the according dimension adjustments into account).  }
  {        However, if this yields no improvement it is advised to analyze }
  {        the integrand in order to determine the integration             }
  {        difficulties. If the position of a local difficulty can be      }
  {        determined (e.g. singularity, discontinuity within the interval)}
  {        one will probably gain from splitting up the interval at this   }
  {        point and calling the integrator on the subranges.              }
  {                                                                        }
  {    = 2 The occurrence of roundoff error is detected, which prevents the}
  {        requested tolerance from being achieved.                        }
  {                                                                        }
  {    = 3 Extremely bad integrand behaviour occurs at some points of the  }
  {        integration interval.                                           }
  {                                                                        }
  {    = 6 The input is invalid, because c = a or c = b, or epsabs <= 0 and}
  {        epsrel < 50*eps_x. result, abserr, neval, last are set to zero. }
  {                                                                        }
  {    = 7 The input is invalid, limit < 0 or limit > QMAXLIM.             }
  {        result, abserr, last are set to zero.                           }
  {                                                                        }
  {    = 8 Dynamic list vectors cannot be allocated.                       }
  {        result, abserr, last are set to zero.                           }
var
  last: integer;
  lsx : word;
  alist, blist, elist, rlist: ^TQXLimArr;
  iord: ^TQILimArr;
begin

  {set error return values}
  result := 0.0;
  abserr := 0.0;
  neval  := 0;

  if limit=0 then limit := DefLimit;
  if (limit<1) or (limit>QMAXLIM) then begin
    ier := 7;
    exit;
  end;

  lsx   := word(limit)*sizeof(extended);
  alist := malloc(lsx);
  blist := malloc(lsx);
  elist := malloc(lsx);
  rlist := malloc(lsx);
  iord  := malloc(limit*sizeof(integer));

  if (alist=nil) or (blist=nil) or (elist=nil) or (rlist=nil) or (iord=nil) then begin
    ier := 8;
  end
  else begin
    qawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ier,
          alist^, blist^, elist^, rlist^, iord^, last);
  end;
  mfree(alist, lsx);
  mfree(blist, lsx);
  mfree(elist, lsx);
  mfree(rlist, lsx);
  mfree(iord,  limit*sizeof(integer));
end;


{---------------------------------------------------------------------}
{---------------------------------------------------------------------}

{ These are my Pascal translations of T. Ooura's DE C/Fortran sources }
{ available from http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html as }
{ http://www.kurims.kyoto-u.ac.jp/~ooura/intde.zip. Original notice:  }
{   Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp) }
{   You may use, copy, modify this code for any purpose and           }
{   without fee. You may distribute this ORIGINAL package.            }


{---------------------------------------------------------------------------}
procedure intde(f: TFuncX; a, b, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }

  {#F}
  { f      - integrand f(x), must be analytic over (a,b)   }
  { a      - lower limit of integration                    }
  { b      - upper limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}
const
  mmax = 256;
  efs  = 0.1;
  hoff = 8.5;
var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ba,ir,iv,h,iback,
  irback,t,ep,em,xw,xa,wg,fa,fb,err,errt,errh,errd: extended;
begin

  result := 0;
  abserr := 0;
  neval  := 0;

  if eps < eps_x/32.0 then begin
    ier := 1;
    exit;
  end
  else begin
    ier := 0;
    if a=b then exit;
  end;

  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ba    := b - a;
  ir    := f(0.5*(a+b));
  ir    := 0.25*ir*ba;
  neval := 1;
  iv    := ir*Pi;
  err   := abs(iv)*epst;
  errh  := 0.0;

  h := 2.0*h0;
  m := 1;
  repeat
    iback  := iv;
    irback := ir;
    t := 0.5*h;
    repeat
      em := exp(t);
      ep := pi_2*em;
      em := pi_2/em;
      repeat
        xw := 1.0/(1.0 + exp(ep-em));
        xa := ba*xw;
        wg := xa*(1.0-xw);
        fa := a+xa;
        fb := b-xa;
        fa := f(fa)*wg;
        fb := f(fb)*wg;
        inc(neval,2);
        ir := ir + (fa+fb);
        iv := iv + (fa+fb)*(ep+em);
        errt:= (abs(fa) + abs(fb))*(ep+em);
        if m=1 then err := err + errt*epst;
        ep := ep*ehp;
        em := em*ehm;
      until (errt <= err) and (xw <= epsh);
      t := t + h;
    until t >= h0;
    if m=1 then begin
      errh := (err/epst)*epsh*h0;
      errd := 1.0 + 2.0*errh;
    end
    else errd := h*(abs(iv - 2.0*iback) + 4.0*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <=errh) or (m > mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else abserr := 0.5*errh*epsh*m/efs;
end;


{---------------------------------------------------------------------------}
procedure intdei(f: TFuncX; a, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { a      - lower limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}
const
  mmax = 256;
  efs  = 0.1;
  hoff = 11.0;
var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ir,iv,h,iback,irback,
  t,ep,em,xp,xm,fp,fm,err,errt,errh,errd: extended;
begin

  result := 0.0;
  abserr := 0.0;
  neval  := 0;

  if eps < eps_x/32.0 then begin
    ier := 1;
    exit;
  end
  else ier := 0;

  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ir    := f(a+1.0);
  neval := 1;
  iv    := ir*pi_2;
  err   := abs(iv)*epst;
  errh  := 0.0;

  h := 2.0*h0;
  m := 1;
  repeat
    iback  := iv;
    irback := ir;
    t := 0.5*h;
    repeat
      em := exp(t);
      ep := pi_4*em;
      em := pi_4/em;
      repeat
        xp  := exp(ep-em);
        xm  := 1.0/xp;
        inc(neval,2);
        fp  := f(a+xp);
        fm  := f(a+xm);
        fp  := fp*xp;
        fm  := fm*xm;
        ir  := ir + (fp+fm);
        iv  := iv + (fp+fm)*(ep+em);
        errt:= (abs(fp) + abs(fm))*(ep+em);
        if m=1 then err := err + errt*epst;
        ep  := ep*ehp;
        em  := em*ehm;
      until (errt <= err) and (xm <= epsh);
      t := t + h;
    until t >= h0;
    if m=1 then begin
      errh := (err/epst)*epsh*h0;
      errd := 1.0 + 2.0*errh;
    end
    else errd := h*(abs(iv - 2.0*iback) + 4.0*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <= errh) or (m >= mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else abserr := 0.5*errh*epsh*m/efs;

end;


{---------------------------------------------------------------------------}
procedure intde_p(f: TFuncXP; p: pointer; a, b, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }

  {#F}
  { f      - integrand f(x), must be analytic over (a,b)   }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { b      - upper limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}
const
  mmax = 256;
  efs  = 0.1;
  hoff = 8.5;
var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ba,ir,iv,h,iback,
  irback,t,ep,em,xw,xa,wg,fa,fb,err,errt,errh,errd: extended;
begin

  result := 0;
  abserr := 0;
  neval  := 0;

  if eps < eps_x/32.0 then begin
    ier   := 1;
    exit;
  end
  else begin
    ier := 0;
    if a=b then exit;
  end;

  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ba    := b - a;
  ir    := f(0.5*(a+b),p);
  ir    := 0.25*ir*ba;
  neval := 1;
  iv    := ir*Pi;
  err   := abs(iv)*epst;
  errh  := 0.0;

  h := 2.0*h0;
  m := 1;
  repeat
    iback  := iv;
    irback := ir;
    t := 0.5*h;
    repeat
      em := exp(t);
      ep := pi_2*em;
      em := pi_2/em;
      repeat
        xw := 1.0/(1.0 + exp(ep-em));
        xa := ba*xw;
        wg := xa*(1.0-xw);
        fa := a+xa;
        fb := b-xa;
        fa := f(fa,p)*wg;
        fb := f(fb,p)*wg;
        inc(neval,2);
        ir := ir + (fa+fb);
        iv := iv + (fa+fb)*(ep+em);
        errt:= (abs(fa) + abs(fb))*(ep+em);
        if m=1 then err := err + errt*epst;
        ep := ep*ehp;
        em := em*ehm;
      until (errt <= err) and (xw <= epsh);
      t := t + h;
    until t >= h0;
    if m=1 then begin
      errh := (err/epst)*epsh*h0;
      errd := 1.0 + 2.0*errh;
    end
    else errd := h*(abs(iv - 2.0*iback) + 4.0*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <=errh) or (m > mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else abserr := 0.5*errh*epsh*m/efs;
end;


{---------------------------------------------------------------------------}
procedure intdei_p(f: TFuncXP; p: pointer; a, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}
const
  mmax = 256;
  efs  = 0.1;
  hoff = 11.0;
var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ir,iv,h,iback,irback,
  t,ep,em,xp,xm,fp,fm,err,errt,errh,errd: extended;
begin

  result := 0.0;
  abserr := 0.0;
  neval  := 0;

  if eps < eps_x/32.0 then begin
    ier := 1;
    exit;
  end
  else ier := 0;

  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ir    := f(a+1.0, p);
  neval := 1;
  iv    := ir*pi_2;
  err   := abs(iv)*epst;
  errh  := 0.0;

  h := 2.0*h0;
  m := 1;
  repeat
    iback  := iv;
    irback := ir;
    t := 0.5*h;
    repeat
      em := exp(t);
      ep := pi_4*em;
      em := pi_4/em;
      repeat
        xp  := exp(ep-em);
        xm  := 1.0/xp;
        inc(neval,2);
        fp  := f(a+xp, p);
        fm  := f(a+xm, p);
        fp  := fp*xp;
        fm  := fm*xm;
        ir  := ir + (fp+fm);
        iv  := iv + (fp+fm)*(ep+em);
        errt:= (abs(fp) + abs(fm))*(ep+em);
        if m=1 then err := err + errt*epst;
        ep  := ep*ehp;
        em  := em*ehm;
      until (errt <= err) and (xm <= epsh);
      t := t + h;
    until t >= h0;
    if m=1 then begin
      errh := (err/epst)*epsh*h0;
      errd := 1.0 + 2.0*errh;
    end
    else errd := h*(abs(iv - 2.0*iback) + 4.0*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <= errh) or (m >= mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else abserr := 0.5*errh*epsh*m/efs;

end;


{---------------------------------------------------------------------------}
procedure intdeo(f: TFuncX; a, omega, eps: extended; var result, abserr: extended;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has an oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { a      - lower limit of integration                    }
  { omega  - frequency of oscillation, i.e. f has the form }
  {          f(x) = g(x)*sin(omega*x + theta) as x -> INF. }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_x/32, eps>1, or omega=0   }
  {          3: max. iterations, result/abserr have values }
  {#F}
const
  mmax  = 256;
  lmax  = 5;
  efs   = 0.1;
  enoff = 0.40;
  pqoff = 2.9;
  ppoff = -0.72;
var
  n,m,l,k: integer;
  epsln,epsh,frq4,per2,pp,pq,ehp,ehm,ir,iv,h,iback,
  irback,t,ep,em,tk,xw,wg,xa,fp,fm,err,errh,tn,errd: extended;
  done: boolean;
begin

  result := 0.0;
  abserr := 0.0;
  neval  := 0;

  if (eps < eps_x/32.0) or (eps > 1.0) or (omega=0.0) then begin
    ier := 1;
    exit;
  end
  else ier := 0;


  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  n     := trunc(enoff*epsln);
  frq4  := abs(omega)/pi_2;
  per2  := Pi/abs(omega);
  pq    := pqoff/epsln;
  pp    := ppoff - ln(pq*pq*frq4);
  ehp   := exp(2.0*pq);
  ehm   := 1.0/ehp;
  xw    := exp(pp - pi_2);
  xa    := a + sqrt(0.5*xw*per2);
  iv    := f(xa);
  neval := 1;
  ir    := iv*xw;
  iv    := 0.5*per2*iv;
  err   := abs(iv);
  errh  := 0.0;

  h := 2.0;
  m := 1;
  repeat
    iback  := iv;
    irback := ir;
    t := 0.5*h;
    repeat
      em := exp(2.0*pq*t);
      ep := pi_4*em;
      em := pi_4/em;
      tk := t;
      repeat
        xw := exp(pp-ep-em);
        wg := sqrt(frq4*xw + tk*tk);
        xa := xw/(tk+wg);
        wg := (pq*xw*(ep-em) + xa)/wg;
        fm := a + xa;
        fp := fm + per2*tk;
        inc(neval,2);
        fm := f(fm);
        fp := f(fp);
        ir := ir + (fp+fm)*xw;
        fm := fm*wg;
        fp := fp*(per2 - wg);
        iv := iv + (fp + fm);
        if m=1 then err := err + (abs(fp) + abs(fm));
        ep := ep*ehp;
        em := em*ehm;
        tk := tk + 1.0;
      until ep >= epsln;

      if m=1 then begin
        errh := err * epsh;
        err  := err * eps;
      end;

      tn := tk;
      while abs(fm) > err do begin
        xw := exp(pp-ep- em);
        xa := 0.5*xw/tk;
        wg := xa*(1.0/tk + 2.0*pq*(ep-em));
        fm := a + xa;
        inc(neval);
        fm := f(fm);
        ir := ir + fm*xw;
        fm := fm*wg;
        iv := iv + fm;
        ep := ep*ehp;
        em := em*ehm;
        tk := tk + 1.0;
      end;

      xa := a + per2*tn;
      inc(neval);
      fm := f(xa);
      em := per2*fm;
      iv := iv + em;
      if (abs(fp) > err) or (abs(em) > err) then begin
        l := 0;
        done := false;
        repeat
          inc(l);
          tn := tn + n;
          em := fm;
          fm := f(a + per2*tn);
          inc(neval);
          xa := fm;
          ep := fm;
          em := em + fm;
          xw := 1;
          wg := 1;
          for k:=1 to pred(n) do begin
            xw := xw*(n+1-k) / k;
            wg := wg + xw;
            fp := f(a + per2*(tn-k));
            inc(neval);
            xa := xa + fp;
            ep := ep + fp*wg;
            em := em + fp*xw;
          end;
          wg := per2*n/(wg*n + xw);
          em := wg*abs(em);
          if (em <= err) or (l >= lmax) then done := true
          else iv := iv + per2*xa;
        until done;
        iv := iv + wg*ep;
        if em>err then err := em;
      end;
      t := t + h;
    until t >= 1.0;
    if m=1 then errd := 1.0 + 2.0*errh
    else errd := h*(abs(iv - 2.0*iback) + pq*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <= errh) or (m >= mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else begin
    ier := 0;
    abserr := 0.5*err*m;
  end;

end;

begin
  DefLimit := 500;
end.
