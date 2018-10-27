function [rr] = hypfun_M_rec0plus_backward(a,b,z,n)
% function [rr] = hypfun_M_rec0plus_backward(a,b,z,n)
% 
% Computes the hypergeometric function \mathbf{M}(a;b-n;z), for n a large
% positive integer, using \mathbf{M}(a;b;z), \mathbf{M}(a;b-1;z) and a
% recurrence relation.
% 
% Copyright John W. Pearson 2014


% Use initial data
a1 = zeros(2,1);
a1(1) = gamfun(b-a)/gamfun(b)*hypergeom(a,b,z);
a1(2) = gamfun(b-a-1)/gamfun(b-1)*hypergeom(a,b-1,z);

% Apply recurrence relation
for j = 2:n
    a1(j+1) = (-(-b+j-z)*a1(j)-z*a1(j-1))/(b-j-a);
end

% Return solution
rr = a1(end)/gamfun(b-a-n);
