function [J] = besfun_taylor(nu,z,tol)
% function [J] = besfun_taylor(nu,z,tol)
% 
% Computes the Bessel function J_{nu}(z), using a Taylor series, up to a
% tolerance tol.
% 
% Copyright John W. Pearson 2014


% Sum series
a1 = 1; b1 = a1;
for j = 1:500
    a1 = (-0.25*z^2)/j/(nu+j)*a1;
    b1 = b1+a1;
    if abs(a1)/abs(b1)<tol
        if imag(nu)>0 || imag(nu)<0 || real(nu)>0
            break
        else
            if j>-real(nu)
                break
            end
        end
    end
    if (j==500) % series has not converged
        [' ' num2str(j) ' terms computed']
        J = (0.5*z)^nu*b1;
        if imag(nu)==0 && imag(z)==0 && real(nu)==fix(real(nu))
            J = real(J);
        end
        return
    end
end

% Compute Bessel function
J = (0.5*z)^nu*b1;

% Special case where nu and z are real
if imag(nu)==0 && imag(z)==0
    if real(z)>0 || (real(z)<0 && real(nu)==fix(real(nu)))
        J = real(J);
    end
end
