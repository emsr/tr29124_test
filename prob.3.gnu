# set terminal png transparent nocrop enhanced size 450,320 font "arial,8" 
# set output 'prob.3.png'
set format x "%.1f" 
set format y "%.1f" 
unset key
set samples 151, 151
set xzeroaxis
set yzeroaxis
set zzeroaxis
set title "arcsin PDF with r = 2.0" 
set xlabel "x ->" 
set xrange [ -3.00000 : 3.00000 ] noreverse nowriteback
set ylabel "probability density ->" 
set yrange [ -0.200000 : 1.20000 ] noreverse nowriteback
isint(x)=(int(x)==x)
Binv(p,q)=exp(lgamma(p+q)-lgamma(p)-lgamma(q))
arcsin(x,r)=r<=0?1/0:abs(x)>r?0.0:invpi/sqrt(r*r-x*x)
beta(x,p,q)=p<=0||q<=0?1/0:x<0||x>1?0.0:Binv(p,q)*x**(p-1.0)*(1.0-x)**(q-1.0)
binom(x,n,p)=p<0.0||p>1.0||n<0||!isint(n)?1/0:  !isint(x)?1/0:x<0||x>n?0.0:exp(lgamma(n+1)-lgamma(n-x+1)-lgamma(x+1)  +x*log(p)+(n-x)*log(1.0-p))
cauchy(x,a,b)=b<=0?1/0:b/(pi*(b*b+(x-a)**2))
chisq(x,k)=k<=0||!isint(k)?1/0:  x<=0?0.0:exp((0.5*k-1.0)*log(x)-0.5*x-lgamma(0.5*k)-k*0.5*log2)
erlang(x,n,lambda)=n<=0||!isint(n)||lambda<=0?1/0:  x<0?0.0:x==0?(n==1?real(lambda):0.0):exp(n*log(lambda)+(n-1.0)*log(x)-lgamma(n)-lambda*x)
extreme(x,mu,alpha)=alpha<=0?1/0:alpha*(exp(-alpha*(x-mu)-exp(-alpha*(x-mu))))
f(x,d1,d2)=d1<=0||!isint(d1)||d2<=0||!isint(d2)?1/0:  Binv(0.5*d1,0.5*d2)*(real(d1)/d2)**(0.5*d1)*x**(0.5*d1-1.0)/(1.0+(real(d1)/d2)*x)**(0.5*(d1+d2))
gmm(x,rho,lambda)=rho<=0||lambda<=0?1/0:  x<0?0.0:x==0?(rho>1?0.0:rho==1?real(lambda):1/0):  exp(rho*log(lambda)+(rho-1.0)*log(x)-lgamma(rho)-lambda*x)
geometric(x,p)=p<=0||p>1?1/0:  !isint(x)?1/0:x<0||p==1?(x==0?1.0:0.0):exp(log(p)+x*log(1.0-p))
halfnormal(x,sigma)=sigma<=0?1/0:x<0?0.0:sqrt2invpi/sigma*exp(-0.5*(x/sigma)**2)
hypgeo(x,N,C,d)=N<0||!isint(N)||C<0||C>N||!isint(C)||d<0||d>N||!isint(d)?1/0:  !isint(x)?1/0:x>d||x>C||x<0||x<d-(N-C)?0.0:exp(lgamma(C+1)-lgamma(C-x+1)-lgamma(x+1)+  lgamma(N-C+1)-lgamma(d-x+1)-lgamma(N-C-d+x+1)+lgamma(N-d+1)+lgamma(d+1)-lgamma(N+1))
laplace(x,mu,b)=b<=0?1/0:0.5/b*exp(-abs(x-mu)/b)
logistic(x,a,lambda)=lambda<=0?1/0:lambda*exp(-lambda*(x-a))/(1.0+exp(-lambda*(x-a)))**2
lognormal(x,mu,sigma)=sigma<=0?1/0:  x<0?0.0:invsqrt2pi/sigma/x*exp(-0.5*((log(x)-mu)/sigma)**2)
maxwell(x,a)=a<=0?1/0:x<0?0.0:fourinvsqrtpi*a**3*x*x*exp(-a*a*x*x)
negbin(x,r,p)=r<=0||!isint(r)||p<=0||p>1?1/0:  !isint(x)?1/0:x<0?0.0:p==1?(x==0?1.0:0.0):exp(lgamma(r+x)-lgamma(r)-lgamma(x+1)+  r*log(p)+x*log(1.0-p))
nexp(x,lambda)=lambda<=0?1/0:x<0?0.0:lambda*exp(-lambda*x)
normal(x,mu,sigma)=sigma<=0?1/0:invsqrt2pi/sigma*exp(-0.5*((x-mu)/sigma)**2)
pareto(x,a,b)=a<=0||b<0||!isint(b)?1/0:x<a?0:real(b)/x*(real(a)/x)**b
poisson(x,mu)=mu<=0?1/0:!isint(x)?1/0:x<0?0.0:exp(x*log(mu)-lgamma(x+1)-mu)
rayleigh(x,lambda)=lambda<=0?1/0:x<0?0.0:lambda*2.0*x*exp(-lambda*x*x)
sine(x,f,a)=a<=0?1/0:  x<0||x>=a?0.0:f==0?1.0/a:2.0/a*sin(f*pi*x/a)**2/(1-sin(twopi*f))
t(x,nu)=nu<0||!isint(nu)?1/0:  Binv(0.5*nu,0.5)/sqrt(nu)*(1+real(x*x)/nu)**(-0.5*(nu+1.0))
triangular(x,m,g)=g<=0?1/0:x<m-g||x>=m+g?0.0:1.0/g-abs(x-m)/real(g*g)
uniform(x,a,b)=x<(a<b?a:b)||x>=(a>b?a:b)?0.0:1.0/abs(b-a)
weibull(x,a,lambda)=a<=0||lambda<=0?1/0:  x<0?0.0:x==0?(a>1?0.0:a==1?real(lambda):1/0):  exp(log(a)+a*log(lambda)+(a-1)*log(x)-(lambda*x)**a)
carcsin(x,r)=r<=0?1/0:x<-r?0.0:x>r?1.0:0.5+invpi*asin(x/r)
cbeta(x,p,q)=ibeta(p,q,x)
cbinom(x,n,p)=p<0.0||p>1.0||n<0||!isint(n)?1/0:  !isint(x)?1/0:x<0?0.0:x>=n?1.0:ibeta(n-x,x+1.0,1.0-p)
ccauchy(x,a,b)=b<=0?1/0:0.5+invpi*atan((x-a)/b)
cchisq(x,k)=k<=0||!isint(k)?1/0:x<0?0.0:igamma(0.5*k,0.5*x)
cerlang(x,n,lambda)=n<=0||!isint(n)||lambda<=0?1/0:x<0?0.0:igamma(n,lambda*x)
cextreme(x,mu,alpha)=alpha<=0?1/0:exp(-exp(-alpha*(x-mu)))
cf(x,d1,d2)=d1<=0||!isint(d1)||d2<=0||!isint(d2)?1/0:1.0-ibeta(0.5*d2,0.5*d1,d2/(d2+d1*x))
cgmm(x,rho,lambda)=rho<=0||lambda<=0?1/0:x<0?0.0:igamma(rho,x*lambda)
cgeometric(x,p)=p<=0||p>1?1/0:  !isint(x)?1/0:x<0||p==0?0.0:p==1?1.0:1.0-exp((x+1)*log(1.0-p))
chalfnormal(x,sigma)=sigma<=0?1/0:x<0?0.0:erf(x/sigma/sqrt2)
chypgeo(x,N,C,d)=N<0||!isint(N)||C<0||C>N||!isint(C)||d<0||d>N||!isint(d)?1/0:  !isint(x)?1/0:x<0||x<d-(N-C)?0.0:x>d||x>C?1.0:hypgeo(x,N,C,d)+chypgeo(x-1,N,C,d)
claplace(x,mu,b)=b<=0?1/0:(x<mu)?0.5*exp((x-mu)/b):1.0-0.5*exp(-(x-mu)/b)
clogistic(x,a,lambda)=lambda<=0?1/0:1.0/(1+exp(-lambda*(x-a)))
clognormal(x,mu,sigma)=sigma<=0?1/0:x<=0?0.0:cnormal(log(x),mu,sigma)
cnormal(x,mu,sigma)=sigma<=0?1/0:0.5+0.5*erf((x-mu)/sigma/sqrt2)
cmaxwell(x,a)=a<=0?1/0:x<0?0.0:igamma(1.5,a*a*x*x)
cnegbin(x,r,p)=r<=0||!isint(r)||p<=0||p>1?1/0:  !isint(x)?1/0:x<0?0.0:ibeta(r,x+1,p)
cnexp(x,lambda)=lambda<=0?1/0:x<0?0.0:1-exp(-lambda*x)
cpareto(x,a,b)=a<=0||b<0||!isint(b)?1/0:x<a?0.0:1.0-(real(a)/x)**b
cpoisson(x,mu)=mu<=0?1/0:!isint(x)?1/0:x<0?0.0:1-igamma(x+1.0,mu)
crayleigh(x,lambda)=lambda<=0?1/0:x<0?0.0:1.0-exp(-lambda*x*x)
csine(x,f,a)=a<=0?1/0:  x<0?0.0:x>a?1.0:f==0?real(x)/a:(real(x)/a-sin(f*twopi*x/a)/(f*twopi))/(1.0-sin(twopi*f)/(twopi*f))
ct(x,nu)=nu<0||!isint(nu)?1/0:0.5+0.5*sgn(x)*(1-ibeta(0.5*nu,0.5,nu/(nu+x*x)))
ctriangular(x,m,g)=g<=0?1/0:  x<m-g?0.0:x>=m+g?1.0:0.5+real(x-m)/g-(x-m)*abs(x-m)/(2.0*g*g)
cuniform(x,a,b)=x<(a<b?a:b)?0.0:x>=(a>b?a:b)?1.0:real(x-a)/(b-a)
cweibull(x,a,lambda)=a<=0||lambda<=0?1/0:x<0?0.0:1.0-exp(-(lambda*x)**a)
gsampfunc(t,n) = t<0?0.5*1/(-t+1.0)**n:1.0-0.5*1/(t+1.0)**n
GPFUN_isint = "isint(x)=(int(x)==x)"
fourinvsqrtpi = 2.25675833419103
invpi = 0.318309886183791
invsqrt2pi = 0.398942280401433
log2 = 0.693147180559945
sqrt2 = 1.4142135623731
sqrt2invpi = 0.797884560802865
twopi = 6.28318530717959
GPFUN_Binv = "Binv(p,q)=exp(lgamma(p+q)-lgamma(p)-lgamma(q))"
GPFUN_arcsin = "arcsin(x,r)=r<=0?1/0:abs(x)>r?0.0:invpi/sqrt(r*r-x*x)"
GPFUN_beta = "beta(x,p,q)=p<=0||q<=0?1/0:x<0||x>1?0.0:Binv(p,q)*x**(p-1.0)*(1.0-x)**(q-1.0)"
GPFUN_binom = "binom(x,n,p)=p<0.0||p>1.0||n<0||!isint(n)?1/0:  !isint(x)?1/0:x<0||x>n?0.0:exp(lgamma(n+1)-lgamma(n-x+1)-lgamma(x+1)  +x*log(p)+(n-x)*log(1.0-p))"
GPFUN_cauchy = "cauchy(x,a,b)=b<=0?1/0:b/(pi*(b*b+(x-a)**2))"
GPFUN_chisq = "chisq(x,k)=k<=0||!isint(k)?1/0:  x<=0?0.0:exp((0.5*k-1.0)*log(x)-0.5*x-lgamma(0.5*k)-k*0.5*log2)"
GPFUN_erlang = "erlang(x,n,lambda)=n<=0||!isint(n)||lambda<=0?1/0:  x<0?0.0:x==0?(n==1?real(lambda):0.0):exp(n*log(lambda)+(n-1.0)*log(x)-lgamma(n)-lambda*x)"
GPFUN_extreme = "extreme(x,mu,alpha)=alpha<=0?1/0:alpha*(exp(-alpha*(x-mu)-exp(-alpha*(x-mu))))"
GPFUN_f = "f(x,d1,d2)=d1<=0||!isint(d1)||d2<=0||!isint(d2)?1/0:  Binv(0.5*d1,0.5*d2)*(real(d1)/d2)**(0.5*d1)*x**(0.5*d1-1.0)/(1.0+(real(d1)/d2)*x)**(0.5*(d1+d2))"
GPFUN_gmm = "gmm(x,rho,lambda)=rho<=0||lambda<=0?1/0:  x<0?0.0:x==0?(rho>1?0.0:rho==1?real(lambda):1/0):  exp(rho*log(lambda)+(rho-1.0)*log(x)-lgamma(rho)-lambda*x)"
GPFUN_geometric = "geometric(x,p)=p<=0||p>1?1/0:  !isint(x)?1/0:x<0||p==1?(x==0?1.0:0.0):exp(log(p)+x*log(1.0-p))"
GPFUN_halfnormal = "halfnormal(x,sigma)=sigma<=0?1/0:x<0?0.0:sqrt2invpi/sigma*exp(-0.5*(x/sigma)**2)"
GPFUN_hypgeo = "hypgeo(x,N,C,d)=N<0||!isint(N)||C<0||C>N||!isint(C)||d<0||d>N||!isint(d)?1/0:  !isint(x)?1/0:x>d||x>C||x<0||x<d-(N-C)?0.0:exp(lgamma(C+1)-lgamma(C-x+1)-lgamma(x+1)+  lgamma(N-C+1)-lgamma(d-x+1)-lgamma(N-C-d+x+1)+lgamma(N-d+1)+lgamma(d+1)-lgamma(N+1))"
GPFUN_laplace = "laplace(x,mu,b)=b<=0?1/0:0.5/b*exp(-abs(x-mu)/b)"
GPFUN_logistic = "logistic(x,a,lambda)=lambda<=0?1/0:lambda*exp(-lambda*(x-a))/(1.0+exp(-lambda*(x-a)))**2"
GPFUN_lognormal = "lognormal(x,mu,sigma)=sigma<=0?1/0:  x<0?0.0:invsqrt2pi/sigma/x*exp(-0.5*((log(x)-mu)/sigma)**2)"
GPFUN_maxwell = "maxwell(x,a)=a<=0?1/0:x<0?0.0:fourinvsqrtpi*a**3*x*x*exp(-a*a*x*x)"
GPFUN_negbin = "negbin(x,r,p)=r<=0||!isint(r)||p<=0||p>1?1/0:  !isint(x)?1/0:x<0?0.0:p==1?(x==0?1.0:0.0):exp(lgamma(r+x)-lgamma(r)-lgamma(x+1)+  r*log(p)+x*log(1.0-p))"
GPFUN_nexp = "nexp(x,lambda)=lambda<=0?1/0:x<0?0.0:lambda*exp(-lambda*x)"
GPFUN_normal = "normal(x,mu,sigma)=sigma<=0?1/0:invsqrt2pi/sigma*exp(-0.5*((x-mu)/sigma)**2)"
GPFUN_pareto = "pareto(x,a,b)=a<=0||b<0||!isint(b)?1/0:x<a?0:real(b)/x*(real(a)/x)**b"
GPFUN_poisson = "poisson(x,mu)=mu<=0?1/0:!isint(x)?1/0:x<0?0.0:exp(x*log(mu)-lgamma(x+1)-mu)"
GPFUN_rayleigh = "rayleigh(x,lambda)=lambda<=0?1/0:x<0?0.0:lambda*2.0*x*exp(-lambda*x*x)"
GPFUN_sine = "sine(x,f,a)=a<=0?1/0:  x<0||x>=a?0.0:f==0?1.0/a:2.0/a*sin(f*pi*x/a)**2/(1-sin(twopi*f))"
GPFUN_t = "t(x,nu)=nu<0||!isint(nu)?1/0:  Binv(0.5*nu,0.5)/sqrt(nu)*(1+real(x*x)/nu)**(-0.5*(nu+1.0))"
GPFUN_triangular = "triangular(x,m,g)=g<=0?1/0:x<m-g||x>=m+g?0.0:1.0/g-abs(x-m)/real(g*g)"
GPFUN_uniform = "uniform(x,a,b)=x<(a<b?a:b)||x>=(a>b?a:b)?0.0:1.0/abs(b-a)"
GPFUN_weibull = "weibull(x,a,lambda)=a<=0||lambda<=0?1/0:  x<0?0.0:x==0?(a>1?0.0:a==1?real(lambda):1/0):  exp(log(a)+a*log(lambda)+(a-1)*log(x)-(lambda*x)**a)"
GPFUN_carcsin = "carcsin(x,r)=r<=0?1/0:x<-r?0.0:x>r?1.0:0.5+invpi*asin(x/r)"
GPFUN_cbeta = "cbeta(x,p,q)=ibeta(p,q,x)"
GPFUN_cbinom = "cbinom(x,n,p)=p<0.0||p>1.0||n<0||!isint(n)?1/0:  !isint(x)?1/0:x<0?0.0:x>=n?1.0:ibeta(n-x,x+1.0,1.0-p)"
GPFUN_ccauchy = "ccauchy(x,a,b)=b<=0?1/0:0.5+invpi*atan((x-a)/b)"
GPFUN_cchisq = "cchisq(x,k)=k<=0||!isint(k)?1/0:x<0?0.0:igamma(0.5*k,0.5*x)"
GPFUN_cerlang = "cerlang(x,n,lambda)=n<=0||!isint(n)||lambda<=0?1/0:x<0?0.0:igamma(n,lambda*x)"
GPFUN_cextreme = "cextreme(x,mu,alpha)=alpha<=0?1/0:exp(-exp(-alpha*(x-mu)))"
GPFUN_cf = "cf(x,d1,d2)=d1<=0||!isint(d1)||d2<=0||!isint(d2)?1/0:1.0-ibeta(0.5*d2,0.5*d1,d2/(d2+d1*x))"
GPFUN_cgmm = "cgmm(x,rho,lambda)=rho<=0||lambda<=0?1/0:x<0?0.0:igamma(rho,x*lambda)"
GPFUN_cgeometric = "cgeometric(x,p)=p<=0||p>1?1/0:  !isint(x)?1/0:x<0||p==0?0.0:p==1?1.0:1.0-exp((x+1)*log(1.0-p))"
GPFUN_chalfnormal = "chalfnormal(x,sigma)=sigma<=0?1/0:x<0?0.0:erf(x/sigma/sqrt2)"
GPFUN_chypgeo = "chypgeo(x,N,C,d)=N<0||!isint(N)||C<0||C>N||!isint(C)||d<0||d>N||!isint(d)?1/0:  !isint(x)?1/0:x<0||x<d-(N-C)?0.0:x>d||x>C?1.0:hypgeo(x,N,C,d)+chypgeo(x-1,N,C,d)"
GPFUN_claplace = "claplace(x,mu,b)=b<=0?1/0:(x<mu)?0.5*exp((x-mu)/b):1.0-0.5*exp(-(x-mu)/b)"
GPFUN_clogistic = "clogistic(x,a,lambda)=lambda<=0?1/0:1.0/(1+exp(-lambda*(x-a)))"
GPFUN_clognormal = "clognormal(x,mu,sigma)=sigma<=0?1/0:x<=0?0.0:cnormal(log(x),mu,sigma)"
GPFUN_cmaxwell = "cmaxwell(x,a)=a<=0?1/0:x<0?0.0:igamma(1.5,a*a*x*x)"
GPFUN_cnegbin = "cnegbin(x,r,p)=r<=0||!isint(r)||p<=0||p>1?1/0:  !isint(x)?1/0:x<0?0.0:ibeta(r,x+1,p)"
GPFUN_cnexp = "cnexp(x,lambda)=lambda<=0?1/0:x<0?0.0:1-exp(-lambda*x)"
GPFUN_cnormal = "cnormal(x,mu,sigma)=sigma<=0?1/0:0.5+0.5*erf((x-mu)/sigma/sqrt2)"
GPFUN_cpareto = "cpareto(x,a,b)=a<=0||b<0||!isint(b)?1/0:x<a?0.0:1.0-(real(a)/x)**b"
GPFUN_cpoisson = "cpoisson(x,mu)=mu<=0?1/0:!isint(x)?1/0:x<0?0.0:1-igamma(x+1.0,mu)"
GPFUN_crayleigh = "crayleigh(x,lambda)=lambda<=0?1/0:x<0?0.0:1.0-exp(-lambda*x*x)"
GPFUN_csine = "csine(x,f,a)=a<=0?1/0:  x<0?0.0:x>a?1.0:f==0?real(x)/a:(real(x)/a-sin(f*twopi*x/a)/(f*twopi))/(1.0-sin(twopi*f)/(twopi*f))"
GPFUN_ct = "ct(x,nu)=nu<0||!isint(nu)?1/0:0.5+0.5*sgn(x)*(1-ibeta(0.5*nu,0.5,nu/(nu+x*x)))"
GPFUN_ctriangular = "ctriangular(x,m,g)=g<=0?1/0:  x<m-g?0.0:x>=m+g?1.0:0.5+real(x-m)/g-(x-m)*abs(x-m)/(2.0*g*g)"
GPFUN_cuniform = "cuniform(x,a,b)=x<(a<b?a:b)?0.0:x>=(a>b?a:b)?1.0:real(x-a)/(b-a)"
GPFUN_cweibull = "cweibull(x,a,lambda)=a<=0||lambda<=0?1/0:x<0?0.0:1.0-exp(-(lambda*x)**a)"
eps = 1e-10
xmin = -3.0
xmax = 3.0
ymin = -5
ymax = 5
GPFUN_gsampfunc = "gsampfunc(t,n) = t<0?0.5*1/(-t+1.0)**n:1.0-0.5*1/(t+1.0)**n"
r = 2.0
mu = 0.0
sigma = 1.41421356237309
plot arcsin(x, r)
