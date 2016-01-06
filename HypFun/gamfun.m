function [g] = gamfun(z)
% function [g] = gamfun(z)
% 
% Master code for computing the Gamma function Gamma(z).
% 
% Uses Godfrey (Lanczos) and Stirling series to compute Gamma function for
% different parameter regimes.
% 
% Reliable up to |z|=171, at which point numerical overflow occurs.
% 
% Copyright John W. Pearson 2014


if real(z)<0.5 % transformation
    g = pi/(sin(pi*z)*gamfun(1-z));
    return
elseif imag(z)==0 && real(z)>0.5 && real(z)==fix(real(z)) % integer case
    g = factorial(z-1); % use built-in MATLAB routine
elseif abs(z)>80
    g = gamfun_stirlingbernoulli(z,eps);
else
    g = gamfun_godfrey(z);
end
