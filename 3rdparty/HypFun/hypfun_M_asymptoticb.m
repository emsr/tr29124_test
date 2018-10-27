function [h] = hypfun_M_asymptoticb(a,b,z,tol)
% function [h] = hypfun_M_asymptoticb(a,b,z,tol)
% 
% Computes the hypergeometric function \mathbf{M}(a;b;z), using an
% asymptotic series (Method (b)), up to a tolerance tol.
% 
% Copyright John W. Pearson 2014


% Initialise quantities
r = zeros(2,1);
A = zeros(2,1);
r(1) = (b-a)*(1-a);
r(2) = (b-a+1)*(2-a)/2;
A(1) = 1+r(1)/z;
A(2) = A(1)+r(1)*r(2)/z^2;

% Compute three-term recurrence repeatedly
for j = 3:500
    % Update r(j) and A(j) in terms of r(j-1), A(j-1) and A(j-2)
    r(j) = (b-a+j-1)*(j-a)/j;
    A(j) = A(j-1)+(A(j-1)-A(j-2))*r(j)/z;
    % Terminate summation if stopping criterion is satisfied
    if abs(A(j)-A(j-1))/abs(A(j-1))<tol && abs(A(j-1)-A(j-2))/abs(A(j-2))<tol
        break
    end
    % No convergence
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Initialise quantities
s = zeros(2,1);
B = zeros(2,1);
s(1) = a*(a-b+1);
s(2) = (a+1)*(a-b+2)/2;
B(1) = 1+s(1)/(-z);
B(2) = B(1)+s(1)*s(2)/z^2;

% Sum series
for j = 3:500
    % Update s(j) and B(j) in terms of s(j), B(j-1) and B(j-2)
    s(j) = (a+j-1)*(a-b+j)/j;
    B(j) = B(j-1)+(B(j-1)-B(j-2))*s(j)/(-z);
    % Stopping criterion
    if abs(B(j)-B(j-1))/abs(B(j-1))<tol && abs(B(j-1)-B(j-2))/abs(B(j-2))<tol
        break
    end
    % No convergence
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

h1 = A(end); h2 = B(end);

% Compute an asymptotic series; which one depends on location of variables
% in the complex plane
if (imag(a)==0 && imag(b)==0 && imag(z)==0 && real(z)>0)
    h3 = real(exp(z)*z^(a-b)/gamfun(a)*h1);
elseif (imag(a)==0 && imag(b)==0 && imag(z)==0 && real(z)<0)
    h3 = real((exp(z)*z^(a-b)/gamfun(a)*h1...
        +exp(-pi*1i*a)*z^(-a)/gamfun(b-a)*h2));
elseif (imag(z)==0)
    h3 = exp(z)*z^(a-b)/gamfun(a)*h1;
elseif (real(z)>0)
    h3 = exp(z)*z^(a-b)/gamfun(a)*h1...
        +exp(pi*1i*a)*z^(-a)/gamfun(b-a)*h2;
elseif (real(z)<0)
    h3 = exp(z)*z^(a-b)/gamfun(a)*h1...
        +exp(-pi*1i*a)*z^(-a)/gamfun(b-a)*h2;
end

% Special case where a, b and z are real
if imag(a)==0 && imag(b)==0 && imag(z)==0
    h = real(h3);
else
    h = h3;
end
