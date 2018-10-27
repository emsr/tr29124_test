function [h] = hypfun_M_asymptotica(a,b,z,tol)
% function [h] = hypfun_M_asymptotica(a,b,z,tol)
% 
% Computes the hypergeometric function \mathbf{M}(a;b;z), using an
% asymptotic series (Method (a)), up to a tolerance tol.
% 
% Copyright John W. Pearson 2014


% Sum series
a1 = 1; b1 = 1;
for j = 1:500
    % Update a1(j) and b1 in terms of previously computed a1(j-1) and b1
    a1 = (b-a+j-1)*(-a+j)/j/z*a1;
    b1 = b1+a1;
    % If stopping criterion is satisfied
    if abs(a1)/abs(b1)<tol
        break
    end
    % No convergence
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Sum series
c1 = 1; d1 = 1;
for k = 1:500
    % Update c1 and d1
    c1 = (a+k-1)*(a-b+k)/k/(-z)*c1;
    d1 = d1+c1;
    % Stopping criterion
    if abs(c1)/abs(d1)<tol
        break
    end
    % Specify if 500 terms computed
    if (k==500)
        [' ' num2str(k) ' terms computed']
        return
    end
end

% Take last terms computed
h1 = b1; h2 = d1;

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
