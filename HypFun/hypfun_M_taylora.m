function [h] = hypfun_M_taylora(a,b,z,tol)
% function [h] = hypfun_M_taylora(a,b,z,tol)
% 
% Computes the hypergeometric function \mathbf{M}(a;b;z), using a Taylor
% series (Method (a)), up to a tolerance tol.
% 
% Copyright John W. Pearson 2014


if (b==fix(b) && b<=0)
    a1 = gamfun(a-b+1)/gamfun(a)/gamfun(-b+2)*z^(-b+1); b1 = a1;
    
    for j = 1:500
        a1 = (a-b+j)/j*z/(j-b+1)*a1;
        b1 = b1+a1;
        % Apply stopping criterion
        if abs(a1)/abs(b1)<tol
            break
        end
        % No convergence
        if (j==500)
            [' ' num2str(j) ' terms computed']
            return
        end
    end
    
    h = b1;
else
    % Initialise relevant quantities
    a1 = 1; b1 = 1;

    for j = 1:500
        % Compute next term
        a1 = (a+j-1)/(b+j-1)*z/j*a1;
        % Update the sum of computed terms
        b1 = b1+a1;
        % Apply stopping criterion
        if abs(a1)/abs(b1)<tol
            break
        end
        % No convergence
        if (j==500)
            [' ' num2str(j) ' terms computed']
            return
        end
    end

    % If convergence achieved, return sum of terms computed
    h = b1/gamfun(b);
end
