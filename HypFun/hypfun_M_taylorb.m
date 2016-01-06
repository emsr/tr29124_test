function [h] = hypfun_M_taylorb(a,b,z,tol)
% function [h] = hypfun_M_taylorb(a,b,z,tol)
% 
% Computes the hypergeometric function \mathbf{M}(a;b;z), using a Taylor
% series (Method (b)), up to a tolerance tol.
% 
% Copyright John W. Pearson 2014


% Initialise vector r, which stores multipliers
r = zeros(2,1);
r(1) = a/b;
r(2) = (a+1)/2/(b+1);

% Initialise vector A, which stores sum of first j+1 terms
A = zeros(2,1);
A(1) = 1+z*r(1);
A(2) = A(1)+z^2*a/b*r(2);

for j=3:500
    % Compute current r_j
    r(j) = (a+j-1)/j/(b+j-1);
    % Compute A_j recursively
    A(j) = A(j-1)+(A(j-1)-A(j-2))*r(j)*z;
    % If stopping criterion is satisfied, terminate summation
    if abs(A(j)-A(j-1))/abs(A(j-1))<tol && abs(A(j-1)-A(j-2))/abs(A(j-2))<tol
        break
    end
    % No convergence
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% If convergence achieved, return sum of terms computed
h = A(end)/gamfun(b);
