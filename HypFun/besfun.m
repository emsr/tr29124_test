function [b] = besfun(nu,z)
% function [b] = besfun(nu,z)
% 
% Master code for computing the Bessel function J_{nu}(z).
% 
% Uses a range of techniques to compute Bessel function for different
% parameter regimes.
% 
% Not complete, but works well for certain regimes of use.
% 
% Copyright John W. Pearson 2014


% Consider different parameter regimes
if real(z)<0
    if imag(nu)==0 && real(nu)<0 && real(nu)==fix(real(nu)) % transformations
        b = besfun(-nu,-z);
    else
        b = (-1)^mu*besfun(nu,-z);
    end
elseif imag(nu)==0 && real(nu)<0 && real(nu)==fix(real(nu)) % transformations
    b = (-1)^nu*besfun(-nu,z);
elseif abs(nu)<30
    if abs(z)<30
        b = besfun_taylor(nu,z,eps);
    else
        b = besfun_asymptotic(nu,z,eps);
    end
elseif abs(nu)>=30
    if abs(z)>=1/30*abs(nu)^2
        b = besfun_asymptotic(nu,z,eps);
    elseif abs(z)<=40 && abs(nu)>=40
        b = besfun_taylor(nu,z,eps);
    else
        ['Another method needed; possibly uniform asymptotics']
        return
    end
end
