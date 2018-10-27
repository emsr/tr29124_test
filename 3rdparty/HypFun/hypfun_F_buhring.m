function [h]=hypfun_F_buhring(a,b,c,z,z0,tol)
% function [h] = hypfun_F_buhring(a,b,c,z,z0,tol)
% 
% Computes the hypergeometric function \mathbf{F}(a,b;c;z), for z near the
% unit disc, using Buhring's expansion formula.
% 
% Copyright John W. Pearson 2014


% Initialise vector for coefficients d_j(a,z0)
d = zeros(2,1);
d(1) = (1+a-b)^(-1)*a*((a+1)*(1-2*z0)+(a+b+1)*z0-c);
d(2) = (2*(2+a-b))^(-1)*(a+1)*(((a+2)*(1-2*z0)+(a+b+1)*z0-c)...
    *(1+a-b)^(-1)*a*((a+1)*(1-2*z0)+(a+b+1)*z0-c)+z0*(1-z0)*a);

% Initialise quantities
a1 = 1+d(1)*(z-z0)^(-1)+d(2)*(z-z0)^(-2);
b1 = a1;
for n = 3:500
    % Update d(n), a1 and b1
    d(n) = (n*(n+a-b))^(-1)*(n+a-1)*...
        (((n+a)*(1-2*z0)+(a+b+1)*z0-c)*d(n-1)+z0*(1-z0)*(n+a-2)*d(n-2));
    a1 = d(n)*(z-z0)^(-n);
    b1 = b1+a1;
    % Convergence test
    if abs(a1)/abs(b1)<tol
        break
    end
    % No convergence
    if (n==500)
        [' ' num2str(n) ' terms computed']
        return
    end
end

% Initialise vector for coefficients e_j(a,z0)
e = zeros(2,1);
e(1) = (1-a+b)^(-1)*b*((b+1)*(1-2*z0)+(a+b+1)*z0-c);
e(2) = (2*(2-a+b))^(-1)*(b+1)*(((b+2)*(1-2*z0)+(a+b+1)*z0-c)...
    *(1-a+b)^(-1)*b*((b+1)*(1-2*z0)+(a+b+1)*z0-c)+z0*(1-z0)*b);

% Initialise quantities
a2 = 1+e(1)*(z-z0)^(-1)+e(2)*(z-z0)^(-2);
b2 = a2;
for n=3:500
    % Update e(n), a2 and b2
    e(n) = (n*(n-a+b))^(-1)*(n+b-1)*...
        (((n+b)*(1-2*z0)+(a+b+1)*z0-c)*e(n-1)+z0*(1-z0)*(n+b-2)*e(n-2));
    a2 = e(n)*(z-z0)^(-n);
    b2 = b2+a2;
    % Convergence test
    if abs(a2)/abs(b2)<tol
        break
    end
    % No convergence
    if (n==500)
        [' ' num2str(n) ' terms computed']
        return
    end
end

% Compute \mathbf{F} using these sums
h = gamfun(b-a)/gamfun(b)/gamfun(c-a)*(z0-z)^(-a)*b1...
    +gamfun(a-b)/gamfun(a)/gamfun(c-b)*(z0-z)^(-b)*b2;
