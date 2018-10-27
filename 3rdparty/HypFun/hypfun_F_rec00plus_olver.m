function [h] = hypfun_F_rec00plus_olver(a,b,c,z,k)
% function [h] = hypfun_F_rec00plus_olver(a,b,c,z,k)
% 
% Computes the hypergeometric function \mathbf{F}(a,b;c+k;z), for k a large
% positive integer, using \mathbf{F}(a,b;c;z) and Olver's algorithm. Note
% that Re(z) must be < 1/2.
% 
% Copyright John W. Pearson 2014


% Set tolerance
tol = eps;

% Apply recurrences to obtain stopping condition
p = [1, -(c+1)*(c-(2*(c+1)-a-b-1)*z)/((c-a+1)*(c-b+1)*z)];
f0 = hypergeom([a,b],c,z);
r = [c*(c+1)*(z-1)*f0/((c-a+1)*(c-b+1)*z), ...
    c*(c+1)^2*(c+2)*(z-1)^2*f0/((c-a+1)*(c-a+2)*(c-b+1)*(c-b+2)*z^2)];
maxpj = 0;
for j = 2:1000
    p(j+1) = -(c+j)*(c+j-1)*(z-1)*p(j-1)/((c-a+j)*(c-b+j)*z)- ...
        (c+j)*(c+j-1-(2*(c+j)-a-b-1)*z)*p(j)/((c-a+j)*(c-b+j)*z);
    r(j+1) = (c+j)*(c+j+1)*(z-1)*r(j)/((c-a+j+1)*(c-b+j+1)*z);
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
h = f(k)/gamfun(c+k);
