function [h] = hypfun_M_rec0plus_olver(a,b,z,k)
% function [h] = hypfun_M_rec0plus_olver(a,b,z,k)
% 
% Computes the hypergeometric function \mathbf{M}(a;b+k;z), for k a large
% positive integer, using \mathbf{M}(a;b;z) and Olver's algorithm.
% 
% Copyright John W. Pearson 2014


% Set tolerance
tol = eps;

% Apply recurrences to obtain stopping condition
p = [1, 1+b/z];
f0 = gamfun(b-a)/gamfun(b)*hypergeom(a,b,z);
r = [(b-a)/z*f0, (b-a)*(b-a+1)/z^2*f0];
maxpj = 0;
for j = 2:1000
    p(j+1) = -(b-a+j-1)*p(j-1)/z+(b+j+z-1)*p(j)/z;
    r(j+1) = (b-a+j)*r(j)/z;
    maxpj = max(maxpj,abs(p(j-1)));
    if abs(r(j)/(p(j)*p(j+1)))*maxpj < tol
        break
    end
end

% Set condition on application of algorithm
N = j;
if k >= N
    ['Method will not work, as k is too large for backward recurrence to be applied']
end

% Apply back substitution
f = zeros(N-1,1);
f(N-1) = r(N-1)/p(N);
for j = N-2:-1:k
    f(j) = (r(j)+p(j)*f(j+1))/p(j+1);
end

% Return solution
h = f(k)/gamfun(b-a+k);
