function [h] = hypfun_M_buchholz(a,b,z,tol)
% function [h] = hypfun_M_buchholz(a,b,z,tol)
% 
% Computes the hypergeometric function \mathbf{M}(a;b;z), using a Buchholz
% series method, up to a tolerance tol.
% 
% Copyright John W. Pearson 2014


if z==real(z) % use MATLAB's Bessel function command for real variable
    % Initialise vector A, which corresponds to entries D_j of series
    A = zeros(3,1);
    A(1) = 1; A(2) = 0; A(3) = b/2;

    % Initialise vectors
    a1 = besselj(b-1,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b-1))...
        +b/2*z^2*besselj(b+1,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b+1));
    b1 = a1;

    for j = 3:500
        % Compute entries of vectors recursively
        A(j+1) = ((j-2+b)*A(j-1)+(2*a-b)*A(j-2))/j;
        a1 = A(j+1)*z^(j)* ...
            besselj(b-1+j,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b-1+j));
        b1 = b1+a1;
        % If stopping criterion is satisfied, terminate computation
        if abs(a1)/abs(b1)<tol
            break
        end
        % No convergence
        if (j==500)
            [' ' num2str(j) ' terms computed']
            return
        end
    end

    % Return computed solution
    h = exp(z/2)*2^(b-1)*b1(end);
else
    % Initialise vector A, which corresponds to entries D_j of series
    A = zeros(3,1);
    A(1) = 1; A(2) = 0; A(3) = b/2;

    % Initialise vectors
    a1 = besfun(b-1,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b-1))...
        +b/2*z^2*besfun(b+1,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b+1));
    b1 = a1;

    for j = 3:500
        % Compute entries of vectors recursively
        A(j+1) = ((j-2+b)*A(j-1)+(2*a-b)*A(j-2))/j;
        a1 = A(j+1)*z^(j)*besfun(b-1+j,sqrt(z*(2*b-4*a)))/(sqrt(z*(2*b-4*a))^(b-1+j));
        b1 = b1+a1;
        % If stopping criterion is satisfied, terminate computation
        if abs(a1)/abs(b1)<tol
            break
        end
        % No convergence
        if (j==500)
            [' ' num2str(j) ' terms computed']
            return
        end
    end

    % Return computed solution
    h = exp(z/2)*2^(b-1)*b1(end);
end
