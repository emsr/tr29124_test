function [rr] = hypfun_M_recplus0_backward(a,b,z,n)
% function [rr] = hypfun_M_recplus0_backward(a,b,z,n)
% 
% Computes the hypergeometric function \mathbf{M}(a-n;b;z), for n a large
% positive integer, using \mathbf{M}(a;b;z), \mathbf{M}(a-1;b;z) and a
% recurrence relation.
%
% Note that this is not the minimal solution of the recurrence relation,
% and hence the solution will not be accurate if n is large.
% 
% Copyright John W. Pearson 2014


% Use initial data
a1 = zeros(2,1);
a1(1) = hypergeom(a,b,z);
a1(2) = hypergeom(a-1,b,z);

% Apply recurrence relation
for j = 2:n
    a1(j+1) = (-(2*(a-j+1)-b+z)*a1(j)+(a-j+1)*a1(j-1))/(b-(a-j+1));
end

% Return solution
rr = a1(end)/gamfun(b);
