function [r] = hypfun_M_recplusplus_miller(a,b,z,n,k)
% function [r] = hypfun_M_rec0plus_miller(a,b,z,n,k)
% 
% Computes the hypergeometric function \mathbf{M}(a+k;b+k;z), for k a large
% positive integer, using \mathbf{M}(a;b;z) and Miller's algorithm, with
% n a number much larger than k.
% 
% Copyright John W. Pearson 2014


% Initialise f and v in Miller's algorithm
f = zeros(k+1,1);
v = zeros(n+1,1);
a1 = zeros(n-1,1);
b1 = zeros(n-1,1);

% Define coefficients in recurrence relation
for i1 = 1:n-1
    a1(i1) = -1/((a+i1)*z);
    b1(i1) = (b+i1-z-1)/((a+i1)*z);
end

% Input minimal solution with k=0
f1 = hypergeom(a,b,z)/gamfun(b);
v(end) = 0; v(end-1) = 1;

% Compute recurrence backwards
for i2 = 2:n
    v(n+1-i2) = -(v(n+3-i2)+b1(n+1-i2)*v(n+2-i2))/a1(n+1-i2);
end

% Apply last line of Miller's algorithm
for i3 = 1:k+1
    f(i3) = f1/v(1)*v(i3)*gamfun(b+i3-1);
end

% Return solution
r = f(end)/gamfun(b+k);
