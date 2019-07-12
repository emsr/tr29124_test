unit sdSDist;

{Double precision special functions: Statistical distributions}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Double precision special functions: Statistical distributions

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

                  [7] Cephes Mathematical Library, Version 2.8
                      http://www.moshier.net/#Cephes or http://www.netlib.org/cephes/
                 [19] Boost C++ Libraries, Release 1.42.0, 2010.
                      http://www.boost.org/
                 [39] R: A Language and Environment for Statistical Computing,
                      Version 2.11.1, http://www.r-project.org/

                  Not used for actual implementation but for background info:

                      H. Rinne, Taschenbuch der Statistik, 4.Auflage, Frankfurt 2008

                      Wikipedia, the free encyclopedia. List of probability distributions
                      with links to many specific probability distributions:
                      http://en.wikipedia.org/wiki/List_of_probability_distributions


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  12.02.13  W.Ehrhardt  Initial BP7 version from AMath.sfSDist
 1.00.01  12.02.13  we          sfd_normstd_cdf, sfd_normstd_pdf

 1.02.00  07.04.13  we          Improved sfd_chi2_pdf for very small x
 1.02.01  13.04.13  we          k changed to longint in sfd_poisson_pmf
 1.02.02  14.04.13  we          sfd_rayleigh_pdf/cdf/inv
 1.02.03  19.04.13  we          sfd_maxwell_pdf/cdf/inv
 1.02.04  20.04.13  we          sfd_evt1_pdf/cdf/inv
 1.02.05  20.04.13  we          sfd_maxwell_cdf/inv with lower igamma functions

 1.05.00  14.08.13  we          Removed code fragment for y>1 in sfd_beta_inv
 1.05.01  14.08.13  we          Check a,b > 0 in sfd_beta_cdf
 1.05.02  15.08.13  we          sfd_kumaraswamy_pdf/cdf/inv
 1.05.03  18.08.13  we          normal/normstd_pdf/cdf with erf_z and erf_p

 1.07.00  19.10.13  we          sfd_moyal_pdf/cdf/inv
 1.07.01  04.11.13  we          improved sfd_t_pdf for small x^2/nu

 1.10.00  04.05.14  we          sfd_zipf_pmf/cdf

 1.11.00  05.06.14  we          sfd_levy_pdf/cdf/inv

 1.12.00  16.06.14  we          code for sdc_ibeta_inv moved to sdGamma2

 1.14.00  07.09.14  we          sfd_logistic_inv uses logit function
 1.14.01  06.10.14  we          sfd_logistic_cdf uses logistic function

 1.16.00  23.12.14  we          sfd_invgamma_pdf/cdf/inv
 1.16.01  26.12.14  we          sfd_ls_cdf/pmf: logarithmic (series) distribution
 1.16.02  03.01.15  we          sfd_wald_cdf/pdf/inv
 1.16.03  05.01.15  we          Error if k<1 in sfd_ls_cdf/pmf

 1.24.00  09.10.16  we          Kolmogorov distribution sfd_ks_cdf/inv

 1.25.00  09.06.17  we          sfd_cauchy_inv uses tanpi

 1.26.00  22.07.17  we          Fix sfd_chi2_pdf for x=0
 1.26.01  22.07.17  we          sfd_chi_pdf/cdf/inv

 1.28.00  06.12.17  we          range checks in sfd_gamma_inv
 1.28.01  06.12.17  we          sfd_nakagami_pdf/cdf/inv
 1.28.02  19.12.17  we          sfd_ks_inv with modified regula falsi

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2009-2017 Wolfgang Ehrhardt

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

(*-------------------------------------------------------------------------
  This Pascal code uses material and ideas from open source and public
  domain libraries, see the file '3rdparty.ama' for the licenses.
---------------------------------------------------------------------------*)


function sfd_beta_pdf(a, b, x: double): double;
  {-Return the probability density function of the beta distribution with}
  { parameters a and b: sfd_beta_pdf = x^(a-1)*(1-x)^(b-1) / beta(a,b)}

function sfd_beta_cdf(a, b, x: double): double;
  {-Return the cumulative beta distribution function, a>0, b>0}

function sfd_beta_inv(a, b, y: double): double;
  {-Return the functional inverse of the beta distribution function. a>0, b>0;}
  { 0 <= y <= 1. Given y the function finds x such that beta_cdf(a, b, x) = y}

function sfd_binomial_cdf(p: double; n, k: longint): double;
  {-Return the cumulative binomial distribution function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}

function sfd_binomial_pmf(p: double; n, k: longint): double;
  {-Return the binomial distribution probability mass function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}

function sfd_cauchy_pdf(a, b, x: double): double;
  {-Return the Cauchy probability density function with location a }
  { and scale b > 0, 1/(Pi*b*(1+((x-a)/b)^2))}

function sfd_cauchy_cdf(a, b, x: double): double;
  {-Return the cumulative Cauchy distribution function with location a}
  { and scale b > 0, = 1/2 + arctan((x-a)/b)/Pi}

function sfd_cauchy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Cauchy distribution function}
  { with location a and scale b > 0}

function sfd_chi_pdf(nu: longint; x: double): double;
  {-Return the probability density function of the chi distribution, nu>0}

function sfd_chi_cdf(nu: longint; x: double): double;
  {-Return the cumulative chi distribution with nu>0 degrees of freedom, x >= 0}

function sfd_chi_inv(nu: longint; p: double): double;
  {-Return the functional inverse of the chi distribution, nu>0, 0 <= p < 1}

function sfd_chi2_pdf(nu: longint; x: double): double;
  {-Return the probability density function of the chi-square distribution, nu>0, x >= 0}

function sfd_chi2_cdf(nu: longint; x: double): double;
  {-Return the cumulative chi-square distribution with nu>0 degrees of freedom, x >= 0}

function sfd_chi2_inv(nu: longint; p: double): double;
  {-Return the functional inverse of the chi-square distribution, nu>0, 0 <= p < 1}

function sfd_evt1_pdf(a, b, x: double): double;
  {-Return the probability density function of the Extreme Value Type I distribution}
  { with location a and scale b > 0, result = exp(-(x-a)/b)/b * exp(-exp(-(x-a)/b)) }

function sfd_evt1_cdf(a, b, x: double): double;
  {-Return the cumulative Extreme Value Type I distribution function}
  { with location a and scale b > 0;  result = exp(-exp(-(x-a)/b)). }

function sfd_evt1_inv(a, b, y: double): double;
  {-Return the functional inverse of the Extreme Value Type I distribution}
  { function with location a and scale b > 0;  result = a - b*ln(ln(-y)). }

function sfd_exp_pdf(a, alpha, x: double): double;
  {-Return the exponential probability density function with location a }
  { and rate alpha > 0, = alpha*exp(-alpha*(x-a)) if x >= a, 0 if x < a.}

function sfd_exp_cdf(a, alpha, x: double): double;
  {-Return the cumulative exponential distribution function with location a}
  { and rate alpha > 0, = 1 - exp(-alpha*(x-a)) if x >= a, 0 if x < a.}

function sfd_exp_inv(a, alpha, y: double): double;
  {-Return the functional inverse of the exponential distribution function with}
  { location a and rate alpha > 0}

function sfd_f_pdf(nu1, nu2: longint; x: double): double;
  {-Return the probability density function of the F distribution; x >= 0, nu1, nu2 > 0}

function sfd_f_cdf(nu1, nu2: longint; x: double): double;
  {-Return the cumulative F distribution function; x >= 0, nu1, nu2 > 0}

function sfd_f_inv(nu1, nu2: longint; y: double): double;
  {-Return the functional inverse of the F distribution, nu1, nu2 > 0, 0 <= y <= 1}

function sfd_gamma_pdf(a, b, x: double): double;
  {-Return the probability density function of a gamma distribution with shape}
  { a>0, scale b>0: sfd_gamma_pdf = x^(a-1)*exp(-x/b)/gamma(a)/b^a, x>0}

function sfd_gamma_cdf(a, b, x: double): double;
  {-Return the cumulative gamma distribution function, shape a>0, scale b>0}

function sfd_gamma_inv(a, b, p: double): double;
  {-Return the functional inverse of the gamma distribution function, shape a>0,}
  { scale b>0, 0 <= p <= 1, i.e. finds x such that gamma_cdf(a, b, x) = p}

function sfd_hypergeo_pmf(n1,n2,n,k: longint): double;
  {-Return the hypergeometric distribution probability mass function; n,n1,n2 >= 0, n <= n1+n2;}
  { i.e. the probability that among n randomly chosen samples from a container}
  { with n1 type1 objects and n2 type2 objects are exactly k type1 objects:}

function sfd_hypergeo_cdf(n1,n2,n,k: longint): double;
  {-Return the cumulative hypergeometric distribution function; n,n1,n2 >= 0, n <= n1+n2}

function sfd_invgamma_pdf(a, b, x: double): double;
  {-Return the probability density function of an inverse gamma distribution with}
  { shape a>0, scale b>0: sfd_invgamma_pdf = (b/x)^a/x*exp(-b/x)/Gamma(a), x >= 0}

function sfd_invgamma_cdf(a, b, x: double): double;
  {-Return the cumulative inverse gamma distribution function, shape a>0, scale b>0:}
  { sfd_invgamma_cdf = Gamma(a,b/x)/Gamma(a) = Q(a,b/x) = igammaq(a,b/x), x >= 0}

function sfd_invgamma_inv(a, b, y: double): double;
  {-Return the functional inverse of the inverse gamma distribution function, shape}
  { a>0, scale b>0, 0 <= y <= 1, i.e. finds x such that invgamma_cdf(a, b, x) = y  }

function sfd_ks_cdf(x: double): double;
  {-Return the limiting form for the cumulative Kolmogorov distribution function}

function sfd_ks_inv(y: double): double;
  {-Return the functional inverse of the Kolmogorov distribution (using bisection)}

function sfd_kumaraswamy_pdf(a, b, x: double): double;
  {-Return the Kumaraswamy probability density function with shape}
  { parameters a,b>0, 0<=x<=1; result = a*b*x^(a-1)*(1-x^a)^(b-1) }

function sfd_kumaraswamy_cdf(a, b, x: double): double;
  {-Return the cumulative Kumaraswamy distribution function with}
  { shape parameters a,b > 0, 0 <= x <= 1; result = 1-(1-x^a)^b}

function sfd_kumaraswamy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Kumaraswamy distribution}
  { with shape parameters a,b > 0; result = [1-(1-y)^(1/b)]^(1/a)}

function sfd_laplace_pdf(a, b, x: double): double;
  {-Return the Laplace probability density function with location a}
  { and scale b > 0, result = exp(-abs(x-a)/b) / (2*b)}

function sfd_laplace_cdf(a, b, x: double): double;
  {-Return the cumulative Laplace distribution function with location a and scale b > 0}

function sfd_laplace_inv(a, b, y: double): double;
  {-Return the functional inverse of the Laplace distribution with location a and scale b > 0}

function sfd_levy_pdf(a, b, x: double): double;
  {-Return the Levy probability density function with}
  { location a and scale parameter b > 0}

function sfd_levy_cdf(a, b, x: double): double;
  {-Return the cumulative Levy distribution function with}
  { location a and scale parameter b > 0}

function sfd_levy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Levy distribution}
  { with location a and scale parameter b > 0}

function sfd_logistic_pdf(a, b, x: double): double;
  {-Return the logistic probability density function with location a}
  { and scale parameter b > 0, exp(-(x-a)/b)/b/(1+exp(-(x-a)/b))^2}

function sfd_logistic_cdf(a, b, x: double): double;
  {-Return the cumulative logistic distribution function with}
  { location a and scale parameter b > 0}

function sfd_logistic_inv(a, b, y: double): double;
  {-Return the functional inverse of the logistic distribution}
  { with location a and scale parameter b > 0}

function sfd_lognormal_pdf(a, b, x: double): double;
  {-Return the log-normal probability density function with}
  { location a and scale parameter b > 0, zero for x <= 0.}

function sfd_lognormal_cdf(a, b, x: double): double;
  {-Return the cumulative log-normal distribution function with}
  { location a and scale parameter b > 0, zero for x <= 0,}

function sfd_lognormal_inv(a, b, y: double): double;
  {-Return the functional inverse of the log-normal distribution}
  { with location a and scale parameter b > 0, 0 < y < 1.}

function sfd_ls_pmf(a: double; k: longint): double;
  {-Return the logarithmic (series) probability mass function}
  { with shape 0 < a < 1, k > 0; result = -a^k/(k*ln(1-a))   }

function sfd_ls_cdf(a: double; k: longint): double;
  {-Return the cumulative logarithmic (series) distribution function with shape 0 < a < 1, k > 0}

function sfd_maxwell_pdf(b, x: double): double;
  {-Return the Maxwell probability density function with scale b > 0, x >= 0}

function sfd_maxwell_cdf(b, x: double): double;
  {-Return the cumulative Maxwell distribution function with scale b > 0, x >= 0}

function sfd_maxwell_inv(b, y: double): double;
  {-Return the functional inverse of the Maxwell distribution with scale b > 0}

function sfd_moyal_pdf(a, b, x: double): double;
  {-Return the Moyal probability density function with}
  { location a and scale parameter b > 0}

function sfd_moyal_cdf(a, b, x: double): double;
  {-Return the cumulative Moyal distribution function with}
  { location a and scale parameter b > 0}

function sfd_moyal_inv(a, b, y: double): double;
  {-Return the functional inverse of the Moyal distribution}
  { with location a and scale parameter b > 0}

function sfd_nakagami_pdf(m, w, x: double): double;
  {-Return the probability density function of the Nakagami distribution with}
  { shape m>0, spread w>0, x>=0: nakagami_pdf = 2x*gamma_pdf(m,w/m,x^2)}

function sfd_nakagami_cdf(m, w, x: double): double;
  {-Return the cumulative Nakagami distribution function, shape m>0, spread w>0}

function sfd_nakagami_inv(m, w, p: double): double;
  {-Return the functional inverse of the Nakagami distribution function, shape m>0,}
  { spread w>0, 0 <= p <= 1, i.e. find x such that nakagami_cdf(m, w, x) = p}

function sfd_negbinom_cdf(p,r: double; k: longint): double;
  {-Return the cumulative negative binomial distribution function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}

function sfd_negbinom_pmf(p,r: double; k: longint): double;
  {-Return the negative binomial distribution probability mass function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}

function sfd_normal_pdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) probability density function with mean mu}
  { and standard deviation sd>0, exp(-0.5*(x-mu)^2/sd^2) / sqrt(2*Pi*sd^2)}

function sfd_normal_cdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) distribution function}
  { with mean mu and standard deviation sd > 0        }

function sfd_normal_inv(mu, sd, y: double): double;
  {-Return the functional inverse of the normal (Gaussian) distribution}
  { with mean mu and standard deviation sd > 0, 0 < y < 1.}

function sfd_normstd_pdf(x: double): double;
  {-Return the std. normal probability density function exp(-x^2/2)/sqrt(2*Pi)}

function sfd_normstd_cdf(x: double): double;
  {-Return the standard normal distribution function}

function sfd_normstd_inv(y: double): double;
  {-Return the inverse standard normal distribution function, 0 < y < 1.}
  { For x=normstd_inv(y) and y from (0,1), normstd_cdf(x) = y}

function sfd_poisson_cdf(mu: double; k: longint): double;
  {-Return the cumulative Poisson distribution function with mean mu >= 0}

function sfd_poisson_pmf(mu: double; k: longint): double;
  {-Return the Poisson distribution probability mass function with mean mu >= 0}

function sfd_pareto_pdf(k, a, x: double): double;
  {-Return the Pareto probability density function with minimum value k > 0}
  { and shape a, x >= a > 0, result = (a/x)*(k/x)^a}

function sfd_pareto_cdf(k, a, x: double): double;
  {-Return the cumulative Pareto distribution function minimum value k > 0}
  { and shape a, x >= a > 0, result = 1-(k/x)^a}

function sfd_pareto_inv(k, a, y: double): double;
  {-Return the functional inverse of the Pareto distribution with minimum}
  { value k > 0 and shape a, x >= a > 0, result = k/(1-x)^(1/a)}

function sfd_rayleigh_pdf(b, x: double): double;
  {-Return the Rayleigh probability density function with}
  { scale b > 0, x >= 0; result = x*exp(-0.5*(x/b)^2)/b^2}

function sfd_rayleigh_cdf(b, x: double): double;
  {-Return the cumulative Rayleigh distribution function with scale b > 0, x >= 0}

function sfd_rayleigh_inv(b, y: double): double;
  {-Return the functional inverse of the Rayleigh distribution with scale b > 0}

function sfd_t_pdf(nu: longint; x: double): double;
  {-Return the probability density function of Student's t distribution, nu>0}

function sfd_t_cdf(nu: longint; t: double): double;
  {-Return the cumulative Student t distribution with nu>0 degrees of freedom}

function sfd_t_inv(nu: longint; p: double): double;
  {-Return the functional inverse of Student's t distribution, nu>0, 0 <= p <= 1}

function sfd_triangular_pdf(a, b, c, x: double): double;
  {-Return the triangular probability density function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}

function sfd_triangular_cdf(a, b, c, x: double): double;
  {-Return the cumulative triangular distribution function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}

function sfd_triangular_inv(a, b, c, y: double): double;
  {-Return the functional inverse of the triangular distribution with}
  { lower limit a, upper limit b, mode c; a<b, a <= c <= b, 0 <= y <= 1}

function sfd_uniform_pdf(a, b, x: double): double;
  {-Return the uniform probability density function on [a,b], a<b}

function sfd_uniform_cdf(a, b, x: double): double;
  {-Return the cumulative uniform distribution function on [a,b], a<b}

function sfd_uniform_inv(a, b, y: double): double;
  {-Return the functional inverse of the uniform distribution on [a,b], a<b}

function sfd_wald_pdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) probability density}
  { function with mean mu > 0, scale b > 0 for x >= 0     }

function sfd_wald_cdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) distribution  }
  { function with mean mu > 0, scale b > 0 for x >= 0}

function sfd_wald_inv(mu, b, y: double): double;
  {-Return the functional inverse of the Wald (inverse Gaussian)}
  { distribution with mean mu > 0, scale b > 0, 0 <= y < 1.     }

function sfd_weibull_pdf(a, b, x: double): double;
  {-Return the Weibull probability density function with shape a > 0}
  { and scale b > 0, result = a*x^(a-1)*exp(-(x/b)^a)/ b^a, x > 0}

function sfd_weibull_cdf(a, b, x: double): double;
  {-Return the cumulative Weibull distribution function with}
  { shape parameter a > 0 and scale parameter b > 0}

function sfd_weibull_inv(a, b, y: double): double;
  {-Return the functional inverse of the Weibull distribution}
  { with shape parameter a > 0 and scale parameter b > 0}

function sfd_zipf_pmf(r: double; k: longint): double;
  {-Return the Zipf distribution probability mass function k^(-(r+1))/zeta(r+1), r>0, k>0}

function sfd_zipf_cdf(r: double; k: longint): double;
  {-Return the cumulative Zipf distribution function H(k,r+1)/zeta(r+1), r>0, k>0}


implementation


uses
  DAMath,
  {$ifopt R+}
    sdBasic,
  {$endif}
  sdErf,
  sdZeta,
  sdGamma,
  sdGamma2;


{---------------------------------------------------------------------------}
function sfd_beta_pdf(a, b, x: double): double;
  {-Return the probability density function of the beta distribution with}
  { parameters a and b: sfd_beta_pdf = x^(a-1)*(1-x)^(b-1) / beta(a,b)}
var
  p,t: double;
begin
  if (a <= 0.0) or (b <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_beta_pdf := NaN_d;
    exit;
  end;
  if (x<0.0) or (x>1.0) then sfd_beta_pdf := 0.0
  else if (x=0.0) or (x=1.0) then begin
    p := power(x, a-1.0);
    p := p*power(1.0-x, b-1.0);
    sfd_beta_pdf := p/sfd_beta(a,b);
  end
  else if b=1.0 then sfd_beta_pdf := a*power(x,a-1.0)
  else begin
    {0 < x < 1}
    t := sfd_ibetaprefix(a,b,x,1-x);
    if t<>0 then begin
      p := x*(1.0-x);
      if p*MaxDouble < t then t := PosInf_d else t := t/p;
    end;
    sfd_beta_pdf := t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_beta_cdf(a, b, x: double): double;
  {-Return the cumulative beta distribution function, a>0, b>0}
begin
  if (a <= 0.0) or (b <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_beta_cdf := NaN_d;
    exit;
  end;
  if x<=0.0 then sfd_beta_cdf := 0.0
  else if x>=1.0 then sfd_beta_cdf := 1.0
  else sfd_beta_cdf := sfd_ibeta(a,b,x);
end;


{---------------------------------------------------------------------------}
function sfd_beta_inv(a, b, y: double): double;
  {-Return the functional inverse of the beta distribution function. a>0, b>0;}
  { 0 <= y <= 1. Given y the function finds x such that beta_cdf(a, b, x) = y}
begin
  sfd_beta_inv := sfd_ibeta_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function sfd_binomial_cdf(p: double; n, k: longint): double;
  {-Return the cumulative binomial distribution function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}
begin
  if (p < 0.0) or (p > 1.0) or (n < 0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_binomial_cdf := NaN_d;
    exit;
  end;
  if k<0 then sfd_binomial_cdf := 0.0
  else if (k>=n) or (p=0.0) then sfd_binomial_cdf := 1.0
  else sfd_binomial_cdf := sfd_ibeta(n-k,k+1,1.0-p)
end;


{---------------------------------------------------------------------------}
function sfd_binomial_pmf(p: double; n, k: longint): double;
  {-Return the binomial distribution probability mass function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}
begin
  if (p < 0.0) or (p > 1.0) or (n < 0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_binomial_pmf := NaN_d;
    exit;
  end;
  if (k<0) or (k>n) then sfd_binomial_pmf := 0.0
  else if n=0 then sfd_binomial_pmf := 1.0
  else if k=0 then sfd_binomial_pmf := power(1.0-p, n)
  else if k=n then sfd_binomial_pmf := power(p, n)
  else sfd_binomial_pmf := sfd_beta_pdf(k+1,n-k+1,p)/(n+1);
end;


{---------------------------------------------------------------------------}
function sfd_cauchy_pdf(a, b, x: double): double;
  {-Return the Cauchy probability density function with location a }
  { and scale b > 0, 1/(Pi*b*(1+((x-a)/b)^2))}
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cauchy_pdf := NaN_d;
    exit;
  end;
  x := sqr((x-a)/b);
  sfd_cauchy_pdf := 1.0/(Pi*b*(1.0+x));
end;


{---------------------------------------------------------------------------}
function sfd_cauchy_cdf(a, b, x: double): double;
  {-Return the cumulative Cauchy distribution function with location a}
  { and scale b > 0, = 1/2 + arctan((x-a)/b)/Pi}
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cauchy_cdf := NaN_d;
    exit;
  end;
  if x=a then sfd_cauchy_cdf := 0.5
  else begin
    {Definition is inaccurate for (x-a)/b < -1, use arctan identity}
    x := x-a;
    if x <= -b then sfd_cauchy_cdf := arctan2(b,-x)/Pi
    else sfd_cauchy_cdf := 0.5 + arctan2(x,b)/Pi;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_cauchy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Cauchy distribution function}
  { with location a and scale b > 0}
begin
  if (b <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cauchy_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_cauchy_inv := PosInf_d
  else if y=0.0 then sfd_cauchy_inv := NegInf_d
  else if y=0.5 then sfd_cauchy_inv := a
  else sfd_cauchy_inv := a - b/tanPi(y);
end;


{---------------------------------------------------------------------------}
function sfd_chi2_pdf(nu: longint; x: double): double;
  {-Return the probability density function of the chi-square distribution, nu>0}
var
  z: double;
const
  xmin  = 2.01e-17; {AMath: 1e-20}
  numax = 39;       {Amath: 448}
begin
  if (nu<1) or (x<0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_chi2_pdf := NaN_d;
    exit;
  end;
  if x<=xmin then begin
    if nu>=numax then sfd_chi2_pdf := 0.0
    else if x=0.0 then begin
      case nu of
          1: sfd_chi2_pdf := posinf_d;
          2: sfd_chi2_pdf := 0.5;
        else sfd_chi2_pdf := 0.0;
      end;
    end
    else begin
      {pdf = (x/2)^(nu/2-1)/(2*GAMMA(nu/2)); with z = sqrt(0.5*x)}
      {pdf = z^(nu-2)/(2*Gamma(nu/2)),  z = sqrt(2*x)/2          }
      z := 0.5*sqrt(2.0*x);
      z := 0.5*power(z,nu-2);
      sfd_chi2_pdf := z/sfd_gamma(0.5*nu);
    end;
  end
  else sfd_chi2_pdf := sfd_igprefix(0.5*nu,0.5*x)/x;
end;


{---------------------------------------------------------------------------}
function sfd_chi2_cdf(nu: longint; x: double): double;
  {-Return the cumulative chi-square distribution with nu>0 degrees of freedom, x >= 0}
begin
  if (nu<1) or (x<0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_chi2_cdf := NaN_d;
    exit;
  end;
  sfd_chi2_cdf := sfd_igammap(0.5*nu, 0.5*x);
end;


{---------------------------------------------------------------------------}
function sfd_chi2_inv(nu: longint; p: double): double;
  {-Return the functional inverse of the chi-square distribution, nu>0, 0 <= p < 1}
begin
  if (nu<1) or (p<0.0) or (p>1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_chi2_inv := NaN_d;
    exit;
  end;
  {will return Inf for p=1}
  sfd_chi2_inv := 2.0*sfd_igammap_inv(0.5*nu, p);
end;


{---------------------------------------------------------------------------}
function sfd_chi_pdf(nu: longint; x: double): double;
  {-Return the probability density function of the chi distribution, nu>0}
const
  p10 = 0.75 + 0.04788456080286535587989; {sqrt(2/Pi)}
begin
  if (x=0.0) and (nu=1) then sfd_chi_pdf := p10
  else sfd_chi_pdf := 2.0*x*sfd_chi2_pdf(nu, x*x);
end;


{---------------------------------------------------------------------------}
function sfd_chi_cdf(nu: longint; x: double): double;
  {-Return the cumulative chi distribution with nu>0 degrees of freedom, x >= 0}
begin
  sfd_chi_cdf := sfd_chi2_cdf(nu, x*x);
end;


{---------------------------------------------------------------------------}
function sfd_chi_inv(nu: longint; p: double): double;
  {-Return the functional inverse of the chi distribution, nu>0, 0 <= p < 1}
begin
  sfd_chi_inv := sqrt(sfd_chi2_inv(nu, p));
end;


{---------------------------------------------------------------------------}
function sfd_evt1_pdf(a, b, x: double): double;
  {-Return the probability density function of the Extreme Value Type I distribution}
  { with location a and scale b > 0, result = exp(-(x-a)/b)/b * exp(-exp(-(x-a)/b)) }
var
  t,z: double;
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_evt1_pdf := NaN_d;
    exit;
  end;
  z := (a-x)/b;
  if z >= ln_MaxDbl then sfd_evt1_pdf := 0.0
  else begin
    t := exp(z);
    z := exp(-t);
    sfd_evt1_pdf := t*z/b;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_evt1_cdf(a, b, x: double): double;
  {-Return the cumulative Extreme Value Type I distribution function}
  { with location a and scale b > 0;  result = exp(-exp(-(x-a)/b)). }
var
  t,z: double;
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_evt1_cdf := NaN_d;
    exit;
  end;
  z := (a-x)/b;
  if z >= ln_MaxDbl then sfd_evt1_cdf := 0.5
  else begin
    t := exp(z);
    sfd_evt1_cdf := exp(-t);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_evt1_inv(a, b, y: double): double;
  {-Return the functional inverse of the Extreme Value Type I distribution}
  { function with location a and scale b > 0;  result = a - b*ln(ln(-y)). }
var
  t: double;
begin
  if (b <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_evt1_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_evt1_inv := PosInf_d
  else if y=0.0 then sfd_evt1_inv := NegInf_d
  else begin
    t := ln(-ln(y));
    sfd_evt1_inv := a - b*t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_exp_pdf(a, alpha, x: double): double;
  {-Return the exponential probability density function with location a }
  { and rate alpha > 0, = alpha*exp(-alpha*(x-a)) if x >= a, 0 if x < a.}
begin
  if alpha <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_exp_pdf := NaN_d;
    exit;
  end;
  if x < a then sfd_exp_pdf := 0.0
  else sfd_exp_pdf := alpha*exp(-alpha*(x-a));
end;


{---------------------------------------------------------------------------}
function sfd_exp_cdf(a, alpha, x: double): double;
  {-Return the cumulative exponential distribution function with location a}
  { and rate alpha > 0, = 1 - exp(-alpha*(x-a)) if x >= a, 0 if x < a.}
begin
  if alpha <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_exp_cdf := NaN_d;
    exit;
  end;
  if x<=a then sfd_exp_cdf := 0.0
  else sfd_exp_cdf := -expm1(-alpha*(x-a));
end;


{---------------------------------------------------------------------------}
function sfd_exp_inv(a, alpha, y: double): double;
  {-Return the functional inverse of the exponential distribution function with}
  { location a and rate alpha > 0}
begin
  if (alpha <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_exp_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_exp_inv := PosInf_d
  else if y=0.0 then sfd_exp_inv := NegInf_d
  else sfd_exp_inv := a - ln1p(-y)/alpha;
end;


{---------------------------------------------------------------------------}
function sfd_f_pdf(nu1, nu2: longint; x: double): double;
  {-Return the probability density function of the F distribution; x >= 0, nu1, nu2 > 0}
var
  a,b,w,z: double;
begin
  if (x<0.0) or (nu1<1) or (nu2<1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_f_pdf := NaN_d;
    exit;
  end;
  a := 0.5*nu1;
  b := 0.5*nu2;
  w := a/b;
  z := (a+b)*ln1p(w*x);
  z := z + sfd_lnbeta(a,b);
  w := a*ln(w) + (a-1.0)*ln(x);
  sfd_f_pdf := exp(w-z);
end;


{---------------------------------------------------------------------------}
function sfd_f_cdf(nu1, nu2: longint; x: double): double;
  {-Return the cumulative F distribution function; x >= 0, nu1, nu2 > 0}
var
  w,a,b: double;
begin
  if (x<0.0) or (nu1<1) or (nu2<1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_f_cdf := NaN_d;
    exit;
  end;
  a := 0.5*nu1;
  b := 0.5*nu2;
  w := a*x;
  w := w/(b+w);
  sfd_f_cdf := sfd_ibeta(a,b,w);
end;


{---------------------------------------------------------------------------}
function sfd_f_inv(nu1, nu2: longint; y: double): double;
  {-Return the functional inverse of the F distribution, nu1, nu2 > 0, 0 <= y <= 1}
var
  w,a,b: double;
begin
  if (y<0.0) or (y>1.0) or (nu1<1) or (nu2<1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_f_inv := NaN_d;
    exit;
  end;
  a := 0.5*nu1;
  b := 0.5*nu2;
  {Set w := 0.5 to get a slightly faster version, but this }
  {is significant only if sfd_beta_inv will be much faster.}
  w := sfd_ibeta(a, b, 0.5);
  if (w > y) or (y>0.999) then begin
    w := sfd_beta_inv(a, b, y);
    sfd_f_inv := (b*w)/(a*(1.0-w));
  end
  else begin
    w := sfd_beta_inv(b, a, 1.0-y);
    sfd_f_inv := (b*(1.0-w))/(a*w);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_gamma_pdf(a, b, x: double): double;
  {-Return the probability density function of a gamma distribution with shape}
  { a>0, scale b>0: sfd_gamma_pdf = x^(a-1)*exp(-x/b)/gamma(a)/b^a, x>0}
begin
  if (a <= 0.0) or (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gamma_pdf := NaN_d;
    exit;
  end;
  if x=0.0 then begin
    if a=1.0 then sfd_gamma_pdf := 1.0/b
    else if a>1.0 then sfd_gamma_pdf := 0
    else sfd_gamma_pdf := PosInf_d
  end
  else sfd_gamma_pdf := sfd_igprefix(a,x/b)/x;
end;


{---------------------------------------------------------------------------}
function sfd_gamma_cdf(a, b, x: double): double;
  {-Return the cumulative gamma distribution function, shape a>0, scale b>0}
begin
  if (a <= 0.0) or (b <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gamma_cdf := NaN_d;
    exit;
  end;
  if x<0.0 then sfd_gamma_cdf := 0.0
  else sfd_gamma_cdf := sfd_igammap(a, x/b);
end;


{---------------------------------------------------------------------------}
function sfd_gamma_inv(a, b, p: double): double;
  {-Return the functional inverse of the gamma distribution function, shape a>0,}
  { scale b>0, 0 <= p <= 1, i.e. finds x such that gamma_cdf(a, b, x) = p}
begin
  if (a <= 0.0) or (b <= 0.0) or (p < 0.0) or (p > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gamma_inv := NaN_d;
    exit;
  end
  else sfd_gamma_inv := sfd_igammap_inv(a, p)*b;
end;


{---------------------------------------------------------------------------}
function bin_raw(p,q: double; n, k: longint): double;
  {-Return raw binomial(n,k)*p^k*q^(n-k) without arg check}
var
  lc,lf: double;

  function bd0(x,np: double): double;
    {-Return "deviance part", see Loader article}
  var
    ej, s, s1, v: double;
    j: integer;
  begin
    if abs(x-np) < 0.1*(x+np) then begin
      v  := (x-np)/(x+np);
      s  := (x-np)*v;
      ej := 2.0*x*v;
      v  := v*v;
      j  := 3;
      repeat
        ej := ej*v;
        s1 := s + ej/j;
        if s=s1 then begin
          bd0 := s1;
          exit;
        end;
        s := s1;
        inc(j,2);
      until j<0;
    end;
    bd0 := x*ln(x/np)+np-x;
  end;

  function stirlerr(k: longint): double;
    {-sfd_lngcorr for integer, avaoid lngamma for small k}
  const
    sterr: array[0..7] of double = (
             0,
             0.81061466795327258219670263595e-1,
             0.41340695955409294093822081405e-1,
             0.27677925684998339148789292755e-1,
             0.20790672103765093111522771755e-1,
             0.16644691189821192163194865355e-1,
             0.13876128823070747998745727025e-1,
             0.11896709945891770095055724095e-1);
  begin
    if k>=8 then stirlerr := sfd_lngcorr(k)
    else stirlerr := sterr[k];
  end;

begin
  {Raw binomial pmf used in hypergeometric pmf based on R function dbinom.c and}
  {Catherine Loader (2000). Fast and Accurate Computation of Binomial Probabilities}
  {http://projects.scipy.org/scipy/raw-attachment/ticket/620/loader2000Fast.pdf}
  if (k<0) or (k>n) then bin_raw := 0.0
  else if n=0 then bin_raw:= 1.0
  else if k=0 then bin_raw:= power(q, n)
  else if k=n then bin_raw:= power(p, n)
  else begin
    lc := stirlerr(n) - stirlerr(k) - stirlerr(n-k) - bd0(k,n*p) - bd0(n-k,n*q);
    lf := ln(TwoPi) + ln(k) + ln1p(-k/n);
    bin_raw := exp(lc - 0.5*lf);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_hypergeo_pmf(n1,n2,n,k: longint): double;
  {-Return the hypergeometric distribution probability mass function; n,n1,n2 >= 0, n <= n1+n2;}
  { i.e. the probability that among n randomly chosen samples from a container}
  { with n1 type1 objects and n2 type2 objects are exactly k type1 objects:}
var
  b1,b2,b3,p,q: double;
begin
  {Result = binomial(n1,k)*binomial(n2,n-k)/binomial(n1+n2,n)}
  if (n1<0) or (n2<0) or (n<0) or (n>n1+n2) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_hypergeo_pmf := NaN_d;
    exit;
  end;
  sfd_hypergeo_pmf := 0.0;
  if (k<0) or (k>n1) or (k>n) then exit;
  if (n2<n) and (k<n-n2) then exit;
  if n1+n2 = 0 then begin
    {here n=n1=n2=0; result is 1 if k=1 and 0 otherwise}
    if k=0 then sfd_hypergeo_pmf := 1.0;
  end
  else begin
    p  := n/(n1+n2);
    q  := (n1+n2-n)/(n1+n2);
    {R trick: use binomial_pmf with p=n/(n1+n2) instead of binomial()}
    b1 := bin_raw(p,q,n1,k);
    b2 := bin_raw(p,q,n2,n-k);
    b3 := bin_raw(p,q,n1+n2,n);
    sfd_hypergeo_pmf := b1*(b2/b3);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_hypergeo_cdf(n1,n2,n,k: longint): double;
  {-Return the cumulative hypergeometric distribution function; n,n1,n2 >= 0, n <= n1+n2}
var
  p,t,r,eps: double;
  i: longint;
begin
  if (n1<0) or (n2<0) or (n<0) or (n>n1+n2) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_hypergeo_cdf := NaN_d;
    exit;
  end;
  if (k >= n1) or (k >= n) then sfd_hypergeo_cdf := 1.0
  else if k<0 then sfd_hypergeo_cdf := 0.0
  else begin
    {Summation of probability mass function using recurrence relations}
    eps := 0.5*eps_d;
    if k >= (n1/(n1+n2))*n then begin
      {use upper tail with ascending recurrence:}
      {f(k+1) = ((n1-k)*(n-k))/((k+1)*(n2-n+k+1)) * f(k)}
      i := succ(k);
      t := sfd_hypergeo_pmf(n1,n2,n,i);
      p := t;
      while (i<n) and (t > eps*p) do begin
        r := ((n1-i)/(i+1)) * ((n-i)/(n2+i+1-n));
        t := t*r;
        p := p+t;
        inc(i);
      end;
      sfd_hypergeo_cdf := 1.0-p;
    end
    else begin
      {use lower tail with descending recurrence:}
      {f(k-1) = (k*(n2-n+k))/((n2-n+k)*(n-k+1)) * f(k)}
      i := k;
      t := sfd_hypergeo_pmf(n1,n2,n,i);
      p := t;
      while (i>0) and (t > eps*p) do begin
        r := (i/(n1-i+1)) * ((n2+i-n)/(n-i+1));
        t := t*r;
        p := p+t;
        dec(i);
      end;
      sfd_hypergeo_cdf := p;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_invgamma_pdf(a, b, x: double): double;
  {-Return the probability density function of an inverse gamma distribution with}
  { shape a>0, scale b>0: sfd_invgamma_pdf = (b/x)^a/x*exp(-b/x)/Gamma(a), x >= 0}
begin
  if (a <= 0.0) or (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_invgamma_pdf := NaN_d;
    exit;
  end;
  if x=0.0 then sfd_invgamma_pdf := 0.0
  else sfd_invgamma_pdf := sfd_igprefix(a,b/x)/x;
end;


{---------------------------------------------------------------------------}
function sfd_invgamma_cdf(a, b, x: double): double;
  {-Return the cumulative inverse gamma distribution function, shape a>0, scale b>0:}
  { sfd_invgamma_cdf = Gamma(a,b/x)/Gamma(a) = Q(a,b/x) = igammaq(a,b/x), x >= 0}
begin
  if (a <= 0.0) or (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_invgamma_cdf := NaN_d;
    exit;
  end;
  if x=0.0 then sfd_invgamma_cdf := 0.0
  else sfd_invgamma_cdf := sfd_igammaq(a, b/x);
end;


{---------------------------------------------------------------------------}
function sfd_invgamma_inv(a, b, y: double): double;
  {-Return the functional inverse of the inverse gamma distribution function, shape}
  { a>0, scale b>0, 0 <= y <= 1, i.e. finds x such that invgamma_cdf(a, b, x) = y  }
begin
  if (a <= 0.0) or (b <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_invgamma_inv := NaN_d;
    exit;
  end;
  if y=0.0 then sfd_invgamma_inv := 0.0
  else if y=1.0 then sfd_invgamma_inv := PosInf_d
  else sfd_invgamma_inv := b/sfd_igammaq_inv(a, y);
end;


{---------------------------------------------------------------------------}
const
  xmax = 4.3662109375;  {8942/2048} {ks_cdf(x)=1}
  xmin = 0.04052734375; {83/2048}   {ks_cdf(x)=0}

{---------------------------------------------------------------------------}
function sfd_ks_cdf(x: double): double;
  {-Return the limiting form for the cumulative Kolmogorov distribution function}
var
  y,z: double;
begin
  {For the algorithms see the AMath manual, formulas from}
  {https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test }
  if x < xmin then sfd_ks_cdf := 0.0
  else if x >= xmax  then sfd_ks_cdf := 1.0
  else if x <= 1.15 then begin
    y := expmx2h(pi_2/x);        {y = exp(-(Pi/x)^2/8}
    z := sqr(sqr(sqr(y)));       {y^8}
    if x > 0.45 then z := z*(1.0 + sqr(z)); {z = y^24 + y^8}
    sfd_ks_cdf := Sqrt_TwoPi*(1.0+z)*y/x;
  end
  else begin
    z := exp(-2.0*sqr(x));
    {sfd_ks_cdf = 1 - 2*(z - z^4 + z^9)}
    if x <= 2.5 then begin
      y := sqr(sqr(z));
      z := z*sqr(y) - y + z;
    end;
    sfd_ks_cdf := 1.0 - 2.0*z;
  end;
end;

{.$define use_damtools}

{---------------------------------------------------------------------------}
function sfd_ks_inv(y: double): double;
  {-Return the functional inverse for Kolmogorov distribution}
var
  x1,x2: double;
{$ifdef use_damtools}
var
  nz, err: integer;
{$else}
var
  x,f,f1,f2,dx,tol: double;
const
  ff = 1.0/1024.0;
{$endif}
begin
  if IsNanD(y) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ks_inv := NaN_d;
    exit;
  end;

  if y=0.0 then begin
    sfd_ks_inv := xmin;
    exit;
  end
  else if y=1.0 then begin
    sfd_ks_inv := xmax;
    exit;
  end;

  {starting values for 5 regions}
  if y >= 0.99 then begin
    x1 := 1.62;
    x2 := xmax;
  end
  else if y >= 0.9 then begin
    x1 := 1.22;
    x2 := 1.63;
  end
  else if y < 1e-5 then begin
    x1 := xmin;
    x2 := 0.31;
  end
  else if y <= 0.1 then begin
    x1 := 0.30;
    x2 := 0.58;
  end
  else begin
    x1 := 0.57;
    x2 := 1.23;
  end;

 {$ifdef use_damtools}

  {If DAMTools can/should be used, the result can be computed with the next}
  {line and the number of sfc_ks_cdf calls will be significantly reduced.  }
  sfd_ks_inv := zbrenty({$ifdef FPC_ProcVar}@{$endif}sfc_ks_cdf, y, x1, x2, 0.0, nz, err);

 {$else}

  tol := 2*eps_d;
  f1  := sfd_ks_cdf(x1)-y;
  f2  := sfd_ks_cdf(x2)-y;
  {Use modified regula falsi loop. This needs only about 40% of the}
  {bisection function calls and is more accurate for small values. }
  {In the loop we always have f1 <= 0, f2 >= 0, f1/(f1-f2) >= 0    }
  repeat
    dx := x2-x1;
    f  := f1/(f1-f2);
    {modification: not too small changes}
    if f < ff then f := ff;
    x  := x1 + dx*f;
    f  := sfd_ks_cdf(x)-y;
    if f < 0.0 then begin
      dx := x1-x;
      x1 := x;
      f1 := f;
    end
    else begin
      dx := x2-x;
      x2 := x;
      f2 := f;
    end;
  until (abs(dx) <= tol*x) or (f=0.0);
  sfd_ks_inv := x;

 {$endif}
end;


{---------------------------------------------------------------------------}
function sfd_kumaraswamy_pdf(a, b, x: double): double;
  {-Return the Kumaraswamy probability density function with shape}
  { parameters a,b>0, 0<=x<=1; result = a*b*x^(a-1)*(1-x^a)^(b-1) }
var
  s,t: double;
begin
  if (a <= 0.0) or (b <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_kumaraswamy_pdf := NaN_d;
    exit;
  end;
  if (x<0.0) or (x>1.0) then sfd_kumaraswamy_pdf := 0.0
  else begin
    s := powm1(x,a);
    t := power(-s, b-1.0);
    s := power(x,a-1.0);
    sfd_kumaraswamy_pdf := a*b*s*t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_kumaraswamy_cdf(a, b, x: double): double;
  {-Return the cumulative Kumaraswamy distribution function with}
  { shape parameters a,b > 0, 0 <= x <= 1; result = 1-(1-x^a)^b}
var
  t: double;
begin
  if (a <= 0.0) or (b <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_kumaraswamy_cdf := NaN_d;
    exit;
  end;
  if x<=0.0 then sfd_kumaraswamy_cdf := 0.0
  else if x>=1.0 then sfd_kumaraswamy_cdf := 1.0
  else begin
    t := powm1(x,a);
    sfd_kumaraswamy_cdf := -powm1(-t,b);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_kumaraswamy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Kumaraswamy distribution}
  { with shape parameters a,b > 0; result = [1-(1-y)^(1/b)]^(1/a)}
var
  t: double;
begin
  if (a <= 0.0) or (b <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_kumaraswamy_inv := NaN_d;
    exit;
  end;
  if (y=0.0) or (y=1.0) then sfd_kumaraswamy_inv := y
  else begin
    t := pow1pm1(-y,1.0/b);
    sfd_kumaraswamy_inv := power(-t,1.0/a);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_laplace_pdf(a, b, x: double): double;
  {-Return the Laplace probability density function with location a}
  { and scale b > 0, result = exp(-abs(x-a)/b) / (2*b)}
var
  y: double;
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_laplace_pdf := NaN_d;
    exit;
  end;
  y := abs(x-a)/b;
  sfd_laplace_pdf := 0.5*exp(-y)/b;
end;


{---------------------------------------------------------------------------}
function sfd_laplace_cdf(a, b, x: double): double;
  {-Return the cumulative Laplace distribution function with location a and scale b > 0}
var
  y: double;
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_laplace_cdf := NaN_d;
    exit;
  end;
  y := (x-a)/b;
  if y>0.0 then sfd_laplace_cdf := 0.5 - 0.5*expm1(-y)
  else sfd_laplace_cdf := 0.5*exp(y);
end;


{---------------------------------------------------------------------------}
function sfd_laplace_inv(a, b, y: double): double;
  {-Return the functional inverse of the Laplace distribution with location a and scale b > 0}
var
  z: double;
begin
  if (b <= 0.0) or (y<0.0) or (y>1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_laplace_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_laplace_inv := PosInf_d
  else if y=0.0 then sfd_laplace_inv := NegInf_d
  else begin
    if y < 0.5 then sfd_laplace_inv := a + b*ln(2.0*y)
    else begin
      z := 2.0*(1.0-y);
      sfd_laplace_inv := a - b*ln(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_levy_pdf(a, b, x: double): double;
  {-Return the Levy probability density function with}
  { location a and scale parameter b > 0}
var
  y,z: double;
const
  zmax = 1491.0; { > -2*ln(succ(0)/2)}
begin
  if b<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_levy_pdf := NaN_d;
    exit;
  end;
  y := x-a;
  if y <= 0.0 then sfd_levy_pdf := 0.0
  else begin
    z := b/y;
    if z>zmax then sfd_levy_pdf := 0.0
    else begin
      y := sqrt(z/TwoPi)/y;
      z := exp(-0.5*z);
      sfd_levy_pdf := y*z;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_levy_cdf(a, b, x: double): double;
  {-Return the cumulative Levy distribution function with}
  { location a and scale parameter b > 0}
var
  y: double;
begin
  if b<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_levy_cdf := NaN_d;
    exit;
  end;
  y := 2.0*(x-a);
  if y<=0.0 then sfd_levy_cdf := 0.0
  else begin
    y := sqrt(b/y);
    sfd_levy_cdf := sfd_erfc(y);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_levy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Levy distribution}
  { with location a and scale parameter b > 0}
var
  z: double;
begin
  if (b<=0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_levy_inv := NaN_d;
    exit;
  end;
  if y=0.0 then sfd_levy_inv := a
  else if y=1.0 then sfd_levy_inv := PosInf_d
  else begin
    z := sfd_erfc_inv(y);
    sfd_levy_inv := a  + 0.5*b/sqr(z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_logistic_pdf(a, b, x: double): double;
  {-Return the logistic probability density function with location a}
  { and scale parameter b > 0, exp(-(x-a)/b)/b/(1+exp(-(x-a)/b))^2}
var
  y,z: double;
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_logistic_pdf := NaN_d;
    exit;
  end;
  z := (x-a)/b;
  if abs(z)>ln_MaxDbl then sfd_logistic_pdf := 0.0
  else begin
    z := exp(z);
    y := z*b;
    if (y>1.0) and (z>MaxDouble/y) then sfd_logistic_pdf := 1.0/y
    else sfd_logistic_pdf := z/b/sqr(1.0+z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_logistic_cdf(a, b, x: double): double;
  {-Return the cumulative logistic distribution function with}
  { location a and scale parameter b > 0}
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_logistic_cdf := NaN_d;
    exit;
  end;
  sfd_logistic_cdf := logistic((x-a)/b);
end;


{---------------------------------------------------------------------------}
function sfd_logistic_inv(a, b, y: double): double;
  {-Return the functional inverse of the logistic distribution}
  { with location a and scale parameter b > 0}
begin
  if (b <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_logistic_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_logistic_inv := PosInf_d
  else if y=0.0 then sfd_logistic_inv := NegInf_d
  else sfd_logistic_inv := a + b*logit(y);
end;


{---------------------------------------------------------------------------}
function sfd_lognormal_pdf(a, b, x: double): double;
  {-Return the log-normal probability density function with}
  { location a and scale parameter b > 0, zero for x <= 0.}
var
  z: double;
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_lognormal_pdf := NaN_d;
    exit;
  end;
  if x <= 0.0 then sfd_lognormal_pdf := 0.0
  else begin
    z := (ln(x) - a)/b;
    x := Sqrt_TwoPi*b*x;
    sfd_lognormal_pdf := exp(-0.5*sqr(z))/x;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_lognormal_cdf(a, b, x: double): double;
  {-Return the cumulative log-normal distribution function with}
  { location a and scale parameter b > 0, zero for x <= 0.}
begin
  if b <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_lognormal_cdf := NaN_d;
    exit;
  end;
  if x<=0 then sfd_lognormal_cdf := 0.0
  else begin
    x := (ln(x)-a)/b;
    sfd_lognormal_cdf := sfd_normstd_cdf(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_lognormal_inv(a, b, y: double): double;
  {-Return the functional inverse of the log-normal distribution}
  { with location a and scale parameter b > 0, 0 < y < 1.}
begin
  if (b <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_lognormal_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_lognormal_inv := PosInf_d
  else if y=0.0 then sfd_lognormal_inv := NegInf_d
  else begin
    y := sfd_normstd_inv(y);
    sfd_lognormal_inv := exp(a + b*y);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ls_pmf(a: double; k: longint): double;
  {-Return the logarithmic (series) probability mass function}
  { with shape 0 < a < 1, k > 0; result = -a^k/(k*ln(1-a))   }
var
  t: double;
begin
  if (a <= 0.0) or (a >= 1.0) or (k<1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ls_pmf := NaN_d;
  end
  else begin
    t := power(a,k)/k;
    sfd_ls_pmf := -t/ln1p(-a);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ls_cdf(a: double; k: longint): double;
  {-Return the cumulative logarithmic (series) distribution function with shape 0 < a < 1, k > 0}
var
  s,t: double;
  i: longint;
begin
  if (a <= 0.0) or (a >= 1.0) or (k<1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ls_cdf := NaN_d;
  end
  else if k<3 then begin
    if k=1 then s := a
    else s := a*(1.0+0.5*a);
    sfd_ls_cdf := -s/ln1p(-a);
  end
  else if (a<=0.5) and (k<=50) then begin
    {Sum pmf}
    t := a;
    s := t;
    i := 2;
    while (i <= k) and (t >= eps_d*s) do begin
      t := t*a;
      s := s + t/i;
      inc(i);
    end;
    sfd_ls_cdf := -s/ln1p(-a);
  end
  else begin
    {Use closed form with Lerch phi}
    t := sfd_lerch(a, 1.0, k);
    s := power(a,k);
    s := s*(t - one_d/k);
    sfd_ls_cdf := 1.0 + s/ln1p(-a);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_maxwell_pdf(b, x: double): double;
  {-Return the Maxwell probability density function with scale b > 0, x >= 0}
var
  y: double;
const
  sqrt2_pi = 0.797884560802865355879892119869;  {sqrt(2/Pi)}
begin
  if (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_maxwell_pdf := NaN_d;
    exit;
  end;
  {pdf = sqrt(2/Pi) * x^2/b^3 * exp(-0.5*(x/b)^2)}
  x := x/b;
  if x >= Sqrt_MaxDbl then sfd_maxwell_pdf := 0.0
  else begin
    y := exp(-0.5*sqr(x));
    y := y*sqr(x)*sqrt2_pi;
    sfd_maxwell_pdf := y/b;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_maxwell_cdf(b, x: double): double;
  {-Return the cumulative Maxwell distribution function with scale b > 0, x >= 0}
var
  y: double;
begin
  if (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_maxwell_cdf := NaN_d;
    exit;
  end;
  {cdf = igammap(3/2,(x/b)^2/2}
  y := x/b;
  if y >= Sqrt_MaxDbl then sfd_maxwell_cdf := 1.0
  else sfd_maxwell_cdf :=sfd_igammap(1.5,0.5*sqr(y));
end;


{---------------------------------------------------------------------------}
function sfd_maxwell_inv(b, y: double): double;
  {-Return the functional inverse of the Maxwell distribution with scale b > 0}
var
  z: double;
begin
  if (b <= 0.0) or (y<0.0) or (y>1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_maxwell_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_maxwell_inv := PosInf_d
  else if y=0.0 then sfd_maxwell_inv := 0.0
  else begin
    {inv = b*sqrt(2*igammap_inv(3/2, y))}
    z := sfd_igammap_inv(1.5, y);
    sfd_maxwell_inv := b*sqrt(2.0*z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_moyal_pdf(a, b, x: double): double;
  {-Return the Moyal probability density function with}
  { location a and scale parameter b > 0}
var
  y,z: double;
begin
  if b<=0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_moyal_pdf := NaN_d;
    exit;
  end;
  z := (x-a)/b;
  if (z<=-14.625) or (z>=2985.0) then sfd_moyal_pdf := 0.0
  else begin
    y := exp(-z);
    y := exp(-0.5*(y+z));
    sfd_moyal_pdf := y/Sqrt_TwoPi/b;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_moyal_cdf(a, b, x: double): double;
  {-Return the cumulative Moyal distribution function with}
  { location a and scale parameter b > 0}
var
  y: double;
begin
  if b<=0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_moyal_cdf := NaN_d;
    exit;
  end;
  y := 0.5*(x-a)/b;
  if y<=-14.625 then sfd_moyal_cdf := 0.0
  else if y>=149.0 then sfd_moyal_cdf := 1.0
  else begin
    y := exp(-y);
    sfd_moyal_cdf := sfd_erfc(y/sqrt2);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_moyal_inv(a, b, y: double): double;
  {-Return the functional inverse of the Moyal distribution}
  { with location a and scale parameter b > 0}
var
  z: double;
begin
  if (b<=0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_moyal_inv := NaN_d;
    exit;
  end;
  if y=0.0 then sfd_moyal_inv := NegInf_d
  else if y=1.0 then sfd_moyal_inv := PosInf_d
  else begin
    z := sfd_erfc_inv(y);
    sfd_moyal_inv := a - b*ln(2.0*sqr(z));
  end;
end;


{---------------------------------------------------------------------------}
function sfd_nakagami_pdf(m, w, x: double): double;
  {-Return the probability density function of the Nakagami distribution with}
  { shape m>0, spread w>0, x>=0: nakagami_pdf = 2x*gamma_pdf(m,w/m,x^2)}
begin
  if (m <= 0.0) or (w <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_nakagami_pdf := NaN_d;
    exit;
  end;
  if x=0.0 then begin
    {$ifdef nakagami_pdf_zero_limit}
      {This would be the limit for x->0}
      if m=0.5 then sfd_nakagami_pdf := 1.0/sqrt(w*Pi_2)
      else if m>0.5 then sfd_nakagami_pdf := 0.0
      else sfd_nakagami_pdf := PosInf_x
    {$else}
      sfd_nakagami_pdf := 0.0;
    {$endif}
  end
  else sfd_nakagami_pdf := sfd_gamma_pdf(m,w/m,x*x)*2.0*x;
end;


{---------------------------------------------------------------------------}
function sfd_nakagami_cdf(m, w, x: double): double;
  {-Return the cumulative Nakagami distribution function, shape m>0, spread w>0}
begin
  if (m <= 0.0) or (w <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_nakagami_cdf := NaN_d;
    exit;
  end;
  if x<0.0 then sfd_nakagami_cdf := 0.0
  else sfd_nakagami_cdf := sfd_igammap(m, sqr(x)*m/w);
end;


{---------------------------------------------------------------------------}
function sfd_nakagami_inv(m, w, p: double): double;
  {-Return the functional inverse of the Nakagami distribution function, shape m>0,}
  { spread w>0, 0 <= p <= 1, i.e. find x such that nakagami_cdf(m, w, x) = p}
begin
  if (p < 0.0) or (p > 1.0) or (m <= 0.0) or (w <= 0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_nakagami_inv := NaN_d;
  end
  else sfd_nakagami_inv := sqrt(sfd_igammap_inv(m, p)*w/m);
end;


{---------------------------------------------------------------------------}
function sfd_negbinom_cdf(p,r: double; k: longint): double;
  {-Return the cumulative negative binomial distribution function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}
begin
  if (p < 0.0) or (p > 1.0) or (r <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_negbinom_cdf := NaN_d;
    exit;
  end;
  if k<0 then sfd_negbinom_cdf := 0.0
  else sfd_negbinom_cdf := sfd_ibeta(r,k+1,p);
end;


{---------------------------------------------------------------------------}
function sfd_negbinom_pmf(p,r: double; k: longint): double;
  {-Return the negative binomial distribution probability mass function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}
var
  t: double;
begin
  {sfd_negbinom_pdf = Gamma(k+r)/(k!*Gamma(r))*p^r*(1-p)^k}
  if (p < 0.0) or (p > 1.0) or (r <= 0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_negbinom_pmf := NaN_d;
    exit;
  end;
  if (k<0) or (p=0.0) then sfd_negbinom_pmf := 0.0
  else begin
    t := sfd_beta_pdf(r,k+1,p);
    sfd_negbinom_pmf := (p/(r+k))*t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_normal_pdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) probability density function with mean mu}
  { and standard deviation sd>0, exp(-0.5*(x-mu)^2/sd^2) / sqrt(2*Pi*sd^2)}
begin
  if sd <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_normal_pdf := NaN_d;
    exit;
  end;
  sfd_normal_pdf := sfd_erf_z((x-mu)/sd)/sd;
end;


{---------------------------------------------------------------------------}
function sfd_normal_cdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) distribution function}
  { with mean mu and standard deviation sd > 0        }
begin
  if sd <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_normal_cdf := NaN_d;
    exit;
  end;
  sfd_normal_cdf := sfd_erf_p((x-mu)/sd);
end;


{---------------------------------------------------------------------------}
function sfd_normal_inv(mu, sd, y: double): double;
  {-Return the functional inverse of the normal (Gaussian) distribution}
  { with mean mu and standard deviation sd > 0, 0 < y < 1.}
begin
  if (sd <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_normal_inv := NaN_d;
    exit;
  end;
  sfd_normal_inv := sd*sfd_normstd_inv(y) + mu;
end;


{---------------------------------------------------------------------------}
function sfd_normstd_pdf(x: double): double;
  {-Return the std. normal probability density function exp(-x^2/2)/sqrt(2*Pi)}
begin
  sfd_normstd_pdf := sfd_erf_z(x);
end;


{---------------------------------------------------------------------------}
function sfd_normstd_cdf(x: double): double;
  {-Return the standard normal distribution function}
begin
  {normstd_cdf = (1 + erf(z))/2 = erfc(z)/2, where z = x/sqrt(2)}
  sfd_normstd_cdf := sfd_erf_p(x);
end;


{---------------------------------------------------------------------------}
function sfd_normstd_inv(y: double): double;
  {-Return the inverse standard normal distribution function, 0 <= y <= 1.}
  { For x=normstd_inv(y) and y from [0,1], normstd_cdf(x) = y}
begin
  if IsNanOrInfD(y) then sfd_normstd_inv := y
  else if y <= 0.0 then sfd_normstd_inv := NegInf_d
  else if y >= 1.0 then sfd_normstd_inv := PosInf_d
  else sfd_normstd_inv := -sfd_erfc_inv(2.0*y)*sqrt2;
end;


{---------------------------------------------------------------------------}
function sfd_pareto_pdf(k, a, x: double): double;
  {-Return the Pareto probability density function with minimum value k > 0}
  { and shape a, x >= a > 0, result = (a/x)*(k/x)^a}
var
  y: double;
begin
  if (k <= 0.0) or (a <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pareto_pdf := NaN_d;
    exit;
  end;
  if x < k then sfd_pareto_pdf := 0.0
  else begin
    y := power(k/x,a);
    sfd_pareto_pdf := y*a/x;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_pareto_cdf(k, a, x: double): double;
  {-Return the cumulative Pareto distribution function minimum value k > 0}
  { and shape a, x >= a > 0, result = 1-(k/x)^a}
begin
  if (k <= 0.0) or (a <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pareto_cdf := NaN_d;
    exit;
  end;
  if x < k then sfd_pareto_cdf := 0.0
  else begin
    sfd_pareto_cdf := -powm1(k/x,a);
    {sfd_pareto_cdf := 1.0 - power(k/x, a);}
    {sfd_pareto_cdf := -expm1(a*ln(k/x));}
  end;
end;


{---------------------------------------------------------------------------}
function sfd_poisson_pmf(mu: double; k: longint): double;
  {-Return the Poisson distribution probability mass function with mean mu >= 0}
var
  t: double;
begin
  if mu<0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_poisson_pmf := NaN_d;
    exit;
  end;
  if k<0.0 then sfd_poisson_pmf := 0.0
  else if mu=0.0 then begin
    if k=0 then sfd_poisson_pmf := 1.0
    else sfd_poisson_pmf := 0.0;
  end
  else begin
    if mu<eps_d then begin
      t := sfd_igprefix(k,mu);
      sfd_poisson_pmf := t/k
    end
    else begin
      t := sfd_igprefix(1.0+k,mu);
      sfd_poisson_pmf := t/mu;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_poisson_cdf(mu: double; k: longint): double;
  {-Return the cumulative Poisson distribution function with mean mu >= 0}
begin
  if mu < 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_poisson_cdf := NaN_d;
    exit;
  end;
  if k<0.0 then sfd_poisson_cdf := 0.0
  else if mu=0.0 then sfd_poisson_cdf := 1.0
  else sfd_poisson_cdf := sfd_igammaq(1.0+k,mu);
end;


{---------------------------------------------------------------------------}
function sfd_pareto_inv(k, a, y: double): double;
  {-Return the functional inverse of the Pareto distribution with minimum}
  { value k > 0 and shape a, x >= a > 0, result = k/(1-x)^(1/a)}
begin
  if (k <= 0.0) or (a <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pareto_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_pareto_inv := PosInf_d
  else if y=0.0 then sfd_pareto_inv := k
  else begin
    sfd_pareto_inv := k/power(1.0-y,1.0/a);
    {sfd_pareto_inv := k*exp(-ln1p(-y)/a);}
  end;
end;


{---------------------------------------------------------------------------}
function sfd_rayleigh_pdf(b, x: double): double;
  {-Return the Rayleigh probability density function with}
  { scale b > 0, x >= 0; result = x*exp(-0.5*(x/b)^2)/b^2}
var
  y: double;
begin
  if (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_rayleigh_pdf := NaN_d;
    exit;
  end;
  x := x/b;
  if x >= Sqrt_MaxDbl then sfd_rayleigh_pdf := 0.0
  else begin
    y := exp(-0.5*sqr(x));
    sfd_rayleigh_pdf := y*x/b;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_rayleigh_cdf(b, x: double): double;
  {-Return the cumulative Rayleigh distribution function with scale b > 0, x >= 0}
var
  y: double;
begin
  if (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_rayleigh_cdf := NaN_d;
    exit;
  end;
  y := x/b;
  if y >= Sqrt_MaxDbl then sfd_rayleigh_cdf := 1.0
  else sfd_rayleigh_cdf := -expm1(-0.5*sqr(y));
end;


{---------------------------------------------------------------------------}
function sfd_rayleigh_inv(b, y: double): double;
  {-Return the functional inverse of the Rayleigh distribution with scale b > 0}
var
  z: double;
begin
  if (b <= 0.0) or (y<0.0) or (y>1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_rayleigh_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_rayleigh_inv := PosInf_d
  else if y=0.0 then sfd_rayleigh_inv := 0.0
  else begin
    z := -2.0*ln1p(-y);
    sfd_rayleigh_inv := b*sqrt(z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_t_pdf(nu: longint; x: double): double;
  {-Return the probability density function of Student's t distribution, nu>0}
var
  b,p,y: double;
begin
  if nu<1 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_t_pdf := NaN_d;
    exit;
  end;
  b := sfd_lnbeta(0.5, 0.5*nu);
  p := x*x/nu;
  {Avoid FPC nonsense}
  y := -0.5; y := y*(nu+1);
  if p<0.0625 then begin
    p := ln1p(p)*y;
    sfd_t_pdf := exp(p-b)/sqrt(nu);
  end
  else begin
    p := power(1.0+p,y);
    sfd_t_pdf := p*exp(-b)/sqrt(nu);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_t_cdf(nu: longint; t: double): double;
  {-Return the cumulative Student t distribution with nu>0 degrees of freedom}
var
  x,z,f,tz,p,xsqk: double;
  j: longint;
begin
  {Ref: Cephes [7], function stdtrl in ldouble/stdtrl.c}
  if nu<1 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_t_cdf := NaN_d;
    exit;
  end;
  if t=0.0 then begin
    sfd_t_cdf := 0.5;
    exit;
  end;
  {Note: Cephes uses only the t < 1.6 condition}
  if (t<-1.6) or (nu>20000) then begin
    if abs(t) >= Sqrt_MaxDbl then z := 0.0 else z := nu/(nu+t*t);
    f := 0.5*sfd_ibeta(0.5*nu, 0.5, z);
    if t<=0.0 then sfd_t_cdf := f else sfd_t_cdf := 1.0-f;
    exit;
  end;
  x := abs(t);
  z := 1.0+x*x/nu;
  if odd(nu) then begin
    {computation for odd nu}
    xsqk := x/sqrt(nu);
    p := arctan(xsqk);
    if nu>1 then begin
      f := 1.0;
      tz := 1.0;
      j := 3;
      while (j <= nu-2) and (tz/f > eps_d) do begin
        tz := tz*((j-1)/(z*j));
        f  := f+tz;
        inc(j,2);
      end;
      p := p + f*xsqk/z;
    end;
    p := p/Pi_2;
  end
  else begin
    f  := 1.0;
    tz := 1.0;
    j  := 2;
    while (j <= nu-2) and (tz/f > eps_d) do begin
      tz := tz*((j-1)/(z*j));
      f  := f+tz;
      inc(j,2);
    end;
    p := f*x/sqrt(z*nu);
  end;
  if t<0.0 then begin
    {note destruction of relative accuracy}
    p := -p;
  end;
  sfd_t_cdf:= 0.5+0.5*p;
end;


{---------------------------------------------------------------------------}
function sfd_t_inv(nu: longint; p: double): double;
  {-Return the functional inverse of Student's t distribution, nu>0, 0 <= p <= 1}
var
  t,rk: double;
  rf: boolean;
begin
  {Ref: Cephes[7], function stdtril in ldouble/stdtrl.c}
  if (nu<1) or (p<0.0) or (p>1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_t_inv := NaN_d;
    exit;
  end;

  if p=0.5 then begin
    sfd_t_inv := 0.0;
    exit;
  end;

  rk := nu;
  rf := p<0.5;

  if (p>0.25) and (p<0.75) then begin
    t := sfd_beta_inv(0.5, 0.5*rk, abs(1.0 - 2.0*p));
    t := sqrt(rk*t/(1.0 - t));
  end
  else begin
    if p>0.5 then p := 1.0-p;
    t := sfd_beta_inv(0.5*rk, 0.5, 2.0*p);
    if t*MaxDouble < rk then t := MaxDouble
    else t := sqrt(rk/t - rk);
  end;
  if rf then sfd_t_inv := -t
  else sfd_t_inv := t;
end;


{---------------------------------------------------------------------------}
function sfd_triangular_pdf(a, b, c, x: double): double;
  {-Return the triangular probability density function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}
var
  t: double;
begin
  if (a >= b) or (c < a) or (c > b) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_triangular_pdf := NaN_d;
    exit;
  end;
  if (x<a) or (x>b) then sfd_triangular_pdf := 0.0
  else begin
    if x=c then sfd_triangular_pdf := 2.0/(b-a)
    else if (x=a) or (x=b) then sfd_triangular_pdf :=0
    else begin
      if x<c then begin
        t := 0.5*(b-a)*(c-a);
        sfd_triangular_pdf := (x-a)/t
      end
      else begin
        t := 0.5*(b-a)*(b-c);
        sfd_triangular_pdf := (b-x)/t
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_triangular_cdf(a, b, c, x: double): double;
  {-Return the cumulative triangular distribution function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}
var
  t: double;
begin
  if (a >= b) or (c < a) or (c > b) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_triangular_cdf := NaN_d;
    exit;
  end;
  if x < a then sfd_triangular_cdf := 0.0
  else if x >= b then sfd_triangular_cdf := 1.0
  else if x=c then sfd_triangular_cdf := (c-a)/(b-a)
  else if x<c then begin
    t := (b-a)*(c-a);
    sfd_triangular_cdf := sqr(x-a)/t;
  end
  else begin
    t := (b-a)*(b-c);
    t := sqr(b-x)/t;
    sfd_triangular_cdf := 1.0 - t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_triangular_inv(a, b, c, y: double): double;
  {-Return the functional inverse of the triangular distribution with}
  { lower limit a, upper limit b, mode c; a<b, a <= c <= b, 0 <= y <= 1}
var
  t: double;
begin
  if (a >= b) or (c < a) or (c > b) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_triangular_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_triangular_inv := b
  else if y=0.0 then sfd_triangular_inv := a
  else begin
    t := (c-a)/(b-a);
    if y=t then sfd_triangular_inv := c
    else if y<t then begin
      t := (b-a)*(c-a)*y;
      sfd_triangular_inv := sqrt(t) + a;
    end
    else begin
      t := (b-a)*(b-c)*(1.0-y);
      sfd_triangular_inv := b-sqrt(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_uniform_pdf(a, b, x: double): double;
  {-Return the uniform probability density function on [a,b], a<b}
begin
  if (a >= b) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_uniform_pdf := NaN_d;
    exit;
  end;
  if (x<a) or (x>b) then sfd_uniform_pdf := 0.0
  else sfd_uniform_pdf := 1.0/(b-a);
end;


{---------------------------------------------------------------------------}
function sfd_uniform_cdf(a, b, x: double): double;
  {-Return the cumulative uniform distribution function on [a,b], a<b}
begin
  if (a >= b) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_uniform_cdf := NaN_d;
    exit;
  end;
  if x < a then sfd_uniform_cdf := 0.0
  else if x >= b then sfd_uniform_cdf := 1.0
  else sfd_uniform_cdf := (x-a)/(b-a);
end;


{---------------------------------------------------------------------------}
function sfd_uniform_inv(a, b, y: double): double;
  {-Return the functional inverse of the uniform distribution on [a,b], a<b}
begin
  if (a >= b) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_uniform_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_uniform_inv := b
  else if y=0.0 then sfd_uniform_inv := a
  else sfd_uniform_inv := a + y*(b-a);
end;


{---------------------------------------------------------------------------}
function sfd_wald_pdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) probability density}
  { function with mean mu > 0, scale b > 0 for x >= 0     }
var
  z,a,t: double;
begin
  if (mu <= 0.0) or (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_wald_pdf := NaN_d;
    exit;
  end;
  if x=0.0 then sfd_wald_pdf := 0.0
  else begin
    z := (x - mu)/mu;
    a := b/x;
    t := sqrt(a)/(x*Sqrt_TwoPi);
    z := (-0.5)*a*sqr(z);
    sfd_wald_pdf := t*exp(z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_wald_cdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) distribution  }
  { function with mean mu > 0, scale b > 0 for x >= 0}
var
  s1,s2,a: double;
begin
  if (mu <= 0.0) or (b <= 0.0) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_wald_cdf := NaN_d;
    exit;
  end;
  if x=0.0 then sfd_wald_cdf := 0.0
  else begin
    s1 := (x - mu)/mu;
    s2 := (x + mu)/mu;
    a  := sqrt(b/x);
    s1 := sfd_erf_p(a*s1);
    s2 := sfd_erf_p(-a*s2);
    a  := exp(2.0*b/mu);
    sfd_wald_cdf := s1 + a*s2;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_wald_inv(mu, b, y: double): double;
  {-Return the functional inverse of the Wald (inverse Gaussian)}
  { distribution with mean mu > 0, scale b > 0, 0 <= y < 1.     }
var
  x,f,d,phi,q: double;
  i,imax: integer;
const
  imax1 = 50;
  imax2 = 500;
begin
  if (mu <= 0.0) or (b <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_wald_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_wald_inv := PosInf_d
  else if y=0.0 then sfd_wald_inv := 0.0
  else begin
    imax := imax1;
    phi  := b/mu;
    if phi > 2.0 then begin
      q := sfd_normstd_inv(y);
      f := q / sqrt(phi) - 0.5/phi;
      x := mu*exp(f);
    end
    else begin
      {This is Boost, the R suppdists code uses qgamma(1.0-y,...)}
      {and this is obviously not appropriate for small y}
      q := sfd_igammaq_inv(0.5, y);
      x := 0.5*b/q;
      if x > 0.5*mu then begin
        q := sfd_igammap_inv(0.5, y);
        f := q / sqrt(phi) - 0.5/phi;
        x := mu*exp(f);
      end;
    end;
    i := 0;
    repeat
      d := sfd_wald_pdf(mu, b, x);
      if d=0.0 then x := 0.0
      else begin
        f := sfd_wald_cdf(mu, b, x);
        q := (y-f)/d;
        x := x + q;
        if x < 0.0 then x := 0.0;
        inc(i);
      end;
      if (x<=0.0) and (imax=imax1) then begin
        {Extreme/invalid values, try restart with x=mode as last option.}
        {This is from Giner & Smyth: "A monotonically convergent Newton }
        {iteration for the quantiles of any unimodal distribution, with }
        {application to the inverse Gaussian distribution". Rev. July 11}
        {2014 of http://www.statsci.org/smyth/pubs/qinvgaussPreprint.pdf}
        {The simple R example code has problems itself, e.g. negative x }
        {for (1,2,0.1) or more than 500 iterations for (1,1,1e-200)     }

        {Accurately compute the mode and increase the iteration limit   }
        f := 1.5/phi;
        f := hypot(1.0,f)+f;
        x := mu/f;
        imax := imax2;
      end;
    until (x=0.0) or (abs(q) <= sqrt_epsh*x) or (i>imax);
    {$ifopt R+}
      if ((x<=0.0) or (i>imax)) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
    {$endif}
    sfd_wald_inv := x;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_weibull_pdf(a, b, x: double): double;
  {-Return the Weibull probability density function with shape a > 0}
  { and scale b > 0, result = a*x^(a-1)*exp(-(x/b)^a)/ b^a, x > 0}
var
  y: double;
begin
  if (b <= 0.0) or (a <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_weibull_pdf := NaN_d;
    exit;
  end;
  if x<0 then sfd_weibull_pdf := 0.0
  else begin
    if x=0.0 then begin
      if a=1.0 then sfd_weibull_pdf := 1.0/b
      else if a>1.0 then sfd_weibull_pdf := 0.0
      else sfd_weibull_pdf := PosInf_d;
    end
    else begin
      y := power(x/b,a);
      y := y*exp(-y);
      sfd_weibull_pdf := y*a/x;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_weibull_cdf(a, b, x: double): double;
  {-Return the cumulative Weibull distribution function with}
  { shape parameter a > 0 and scale parameter b > 0}
var
  y: double;
begin
  if (b <= 0.0) or (a <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_weibull_cdf := NaN_d;
    exit;
  end;
  if x<=0.0 then sfd_weibull_cdf := 0.0
  else begin
    y := power(x/b, a);
    sfd_weibull_cdf := -expm1(-y);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_weibull_inv(a, b, y: double): double;
  {-Return the functional inverse of the Weibull distribution}
  { with shape parameter a > 0 and scale parameter b > 0}
var
  z: double;
begin
  if (b <= 0.0) or (a <= 0.0) or (y < 0.0) or (y > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_weibull_inv := NaN_d;
    exit;
  end;
  if y=1.0 then sfd_weibull_inv := PosInf_d
  else if y=0.0 then sfd_weibull_inv := 0
  else begin
    z := -ln1p(-y);
    sfd_weibull_inv := b*power(z,1.0/a);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_zipf_pmf(r: double; k: longint): double;
  {-Return the Zipf distribution probability mass function k^(-(r+1))/zeta(r+1), r>0, k>0}
var
  p,z: double;
begin
  if r <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_zipf_pmf := NaN_d;
  end
  else if k<1 then sfd_zipf_pmf := 0.0
  else begin
    z := sfd_zeta1p(r);
    p := power(k,-(r+1.0));
    sfd_zipf_pmf := p/z;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_zipf_cdf(r: double; k: longint): double;
  {-Return the cumulative Zipf distribution function H(k,r+1)/zeta(r+1), r>0, k>0}
var
  h,z: double;
begin
  if r <= 0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_zipf_cdf := NaN_d;
  end
  else if k<1 then sfd_zipf_cdf := 0.0
  else begin
    {sfd_zipf_cdf := sfd_harmonic2(k,r+1.0)/sfd_zeta(r+1.0);}
    z := sfd_zeta1p(r);
    {Avoid double computation of zeta(r+1) using property of H}
    {H(k,r+1) = zeta(r+1) - zetah(r+1,k+1) and therefore      }
    {sfd_zipf_cdf = 1 - zetah(r+1, k+1) / zeta(r+1)           }
    h := sfd_zetah(r+1.0, k+1);
    sfd_zipf_cdf := 1.0 - h/z;
  end;
end;

end.
