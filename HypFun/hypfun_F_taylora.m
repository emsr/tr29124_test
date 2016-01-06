function [h] = hypfun_F_taylora(a,b,c,z,tol)
% function [h] = hypfun_F_taylora(a,b,c,z,tol)
% 
% Computes the hypergeometric function \mathbf{F}(a,b;c;z), using a Taylor
% series (Method (a)), up to a tolerance tol.
% 
% Copyright John W. Pearson 2014


if (c==fix(c) && c<=0)
    a1 = gamfun(a-c+1)/gamfun(a)*gamfun(b-c+1)/gamfun(b)/gamfun(-c+2)*z^(-c+1);
    b1 = a1;
    
    for j = 1:500
        a1 = (a-c+j)*(b-c+j)/j*z/(j-c+1)*a1;
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

    for j=1:500
        % Compute next term
        a1 = (a+j-1)*(b+j-1)/(c+j-1)*z/j*a1;
        b1 = b1+a1;
        % Update the sum of computed terms
        if abs(a1)/abs(b1)<tol
            break
        end
        % Apply stopping criterion
        if (j==500)
            [' ' num2str(j) ' terms computed']
            return
        end
    end

    % If convergence achieved, return sum of terms computed
    h = b1/gamfun(c);
end 
