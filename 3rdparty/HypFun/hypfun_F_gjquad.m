function [h] = hypfun_F_gjquad(a,b,c,z,N)
% function [h] = hypfun_F_gjquad(a,b,c,z,N)
% 
% Computes the hypergeometric function \mathbf{F}(a,b;c;z), using
% Gauss-Jacobi quadrature with N points.
% 
% Copyright John W. Pearson 2014


% Compute Gauss-Jacobi quadrature nodes and weights
[x,w] = gjquad_nw(c-b-1,b-1,N);

% Use these to compute \mathbf{F} using its integral representation
h = 1/gamfun(b)/gamfun(c-b)/(2^(c-1))*sum(w.*(1-z/2*(x+1)).^(-a));
