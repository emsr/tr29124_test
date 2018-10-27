function [rr] = hypfun_F_rec00plus_backward(a,b,c,z,n)
% function [rr] = hypfun_F_rec00plus_backward(a,b,c,z,n)
% 
% Computes the hypergeometric function \mathbf{F}(a,b;c-n;z), for n a large
% positive integer, using \mathbf{F}(a,b;c;z), \mathbf{F}(a,b;c-1;z) and a
% recurrence relation.
% 
% Copyright John W. Pearson 2014


% Use initial data
a1 = zeros(2,1);
a1(1) = hypergeom([a,b],c,z);
a1(2) = hypergeom([a,b],c-1,z);

% Apply recurrence relation
for j = 2:n
    a1(j+1)=(-(c+1-j)*(c-j-(2*c+1-2*j-a-b)*z)*a1(j)- ...
        (c-a+1-j)*(c-b+1-j)*z*a1(j-1))/(c+1-j)/(c-j)/(z-1);
end

% Return solution
rr = a1(end)/gamfun(c-n);
