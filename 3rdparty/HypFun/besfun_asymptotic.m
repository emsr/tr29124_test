function [J] = besfun_asymptotic(nu,z,tol)
% function [b] = besfun_asymptotic(nu,z,tol)
% 
% Computes the Bessel function J_{nu}(z), using an asymptotic series, up to
% a tolerance tol.
% 
% Copyright John W. Pearson 2014


% Compute parameters
xi = z-(0.5*nu+0.25)*pi;
mu = 4*nu^2;

% Sum first series
a1 = 1; P = a1;
for j = 1:500
    a1 = -(mu-(4*j-3)^2)*(mu-(4*j-1)^2)/(2*j-1)/(2*j)/(8*z)^2*a1;
    P = P+a1;
    if abs(a1)/abs(P)<tol && j>0.5*nu-0.25
        break
    end
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Sum second series
b1 = (mu-1)/(8*z); Q = b1;
for j = 1:500
    b1 = -(mu-(4*j-1)^2)*(mu-(4*j+1)^2)/(2*j)/(2*j+1)/(8*z)^2*b1;
    Q = Q+b1;
    if abs(b1)/abs(Q)<tol && j>0.5*nu-0.75
        break
    end
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Compute Bessel function
J = sqrt(2/(pi*z))*(P*cos(xi)-Q*sin(xi));

% Special case where nu and z are real
if imag(nu)==0 && imag(z)==0 && real(nu)==fix(real(nu))
    J = real(J);
end
