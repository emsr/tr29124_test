function [h] = hypfun_M_singlefraction(a,b,z,tol)
% function [h] = hypfun_M_singlefraction(a,b,z,tol)
% 
% Computes the hypergeometric function \mathbf{M}(a;b;z), using a Taylor
% series (Method (c)), the single fraction method, up to a tolerance tol.
% 
% Copyright John W. Pearson 2014


% Initialise vectors
a1 = [0,b]; b1 = [1,a*z]; c1 = [1,b]; d1 = [1,(b+a*z)/b];

for j = 3:500
    % Compute next values of a1(j),b1(j),c1(j) recursively
    a1(j) = (a1(j-1)+b1(j-1))*(j-1)*(b+j-2);
    b1(j) = b1(j-1)*(a+j-2)*z;
    c1(j) = c1(j-1)*(j-1)*(b+j-2);
    % Stop if any of these are infinitely large
    if (a1(j)==Inf) || (b1(j)==Inf) || (c1(j)==Inf)
        break
    end
    % Compute next term in d1(j)
    d1(j) = (a1(j)+b1(j))/c1(j);
    % Apply stopping criterion
    if abs(d1(j)-d1(j-1))/abs(d1(j-1))<tol && abs(d1(j-1)-d1(j-2))/abs(d1(j-2))<tol
        break
    end
    % No convergence
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% If convergence achieved, return sum of terms computed
h = d1(end)/gamfun(b);
