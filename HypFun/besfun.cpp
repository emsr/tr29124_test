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
if (std::real(z) < 0)
{
    if (std::imag(nu) == 0 && std::real(nu) < 0 && std::real(nu) == fix(real(nu))) // transformations
        b = besfun(-nu, -z);
    else
        b = (-1)^mu * besfun(nu, -z);
}
else if (std::imag(nu) == 0 && std::real(nu)<0 && std::real(nu) == fix(std::real(nu))) // transformations
    b = (-1)^nu * besfun(-nu, z);
else if (std::abs(nu) < 30)
{
    if (std::abs(z) < 30)
        b = besfun_taylor(nu, z, eps);
    else
        b = besfun_asymptotic(nu, z, eps);
}
else if (std::abs(nu) >= 30)
{
    if (std::abs(z) >= 1/30 * std::abs(nu)^2)
        b = besfun_asymptotic(nu, z, eps);
    else if (std::abs(z) <= 40 && std::abs(nu) >= 40)
        b = besfun_taylor(nu, z, eps);
    else//['Another method needed; possibly uniform asymptotics']
        return;
}
