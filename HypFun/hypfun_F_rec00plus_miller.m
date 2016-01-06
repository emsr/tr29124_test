function [r] = hypfun_F_rec00plus_miller(a,b,c,z,n,k)
% function [r] = hypfun_F_rec00plus_miller(a,b,c,z,n,k)
% 
% Computes the hypergeometric function \mathbf{F}(a,b;c+k;z), for k a large
% positive integer, using \mathbf{F}(a,b;c;z) and Miller's algorithm, with
% n a number much larger than k. Note that Re(z) must be < 1/2.
% 
% Copyright John W. Pearson 2014


% Initialise f and v in Miller's algorithm
f = zeros(k+1,1);
v = zeros(n+1,1);
a1 = zeros(n-1,1);
b1 = zeros(n-1,1);

% Define coefficients in recurrence relation
for i1 = 1:n-1
    a1(i1) = (c+i1)*(c+i1-1)*(z-1)/(c+i1-a)/(c+i1-b)/z;
    b1(i1) = (c+i1)*(c+i1-1-(2*(c+i1)-a-b-1)*z)/(c+i1-a)/(c+i1-b)/z;
end

% Input minimal solution with k=0
f1 = hypergeom([a,b],c,z);
v(end) = 0; v(end-1) = 1;

% Compute recurrence backwards
for i2 = 2:n
    v(n+1-i2) = -(v(n+3-i2)+b1(n+1-i2)*v(n+2-i2))/a1(n+1-i2);
end

% Apply last line of Miller's algorithm
for i3 = 1:k+1
    f(i3) = f1/v(1)*v(i3);
end

% Return solution
r = f(end)/gamfun(c+k);
