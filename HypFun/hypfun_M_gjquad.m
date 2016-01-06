function [h] = hypfun_M_gjquad(a,b,z,N)
% function [h] = hypfun_M_gjquad(a,b,z,tol)
% 
% Computes the hypergeometric function \mathbf{M}(a;b;z), using
% Gauss-Jacobi quadrature with N points.
% 
% Copyright John W. Pearson 2014


% Compute Gauss-Jacobi quadrature nodes and weights
[x,w] = gjquad_nw(b-a-1,a-1,N);

% Use these to compute \mathbf{M} using its integral representation
h = 1/gamfun(a)/gamfun(b-a)/(2^(b-1))*exp(0.5*z)*sum(w.*(exp(0.5*z*x)));
