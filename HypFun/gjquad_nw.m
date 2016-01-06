function [nodes,weights] = gjquad_nw(alpha,beta,N)
% function [nodes,weights] = gjquad_nw(alpha,beta,N)
% 
% Computes Gauss-Jacobi nodes and weights using Golub-Welsch method.
% 
% This code uses the notation of Golub and Welsch.
% 
% Copyright John W. Pearson 2014


% Apply condition for quadrature rule to be valid
if real(alpha)<=-1 || real(beta)<=-1,
    ['Re(alpha) and Re(beta) must both be >= -1']
    return
end

% Initialize vectors
a = zeros(N,1);
b = zeros(N-1,1);

% Compute crucial quantities
ab2 = 2+alpha+beta;
a(1) = (beta-alpha)/ab2;
b(1) = sqrt((4*(1+alpha)*(1+beta))/((ab2+1)*ab2^2));

% Compute entries of vectors a and b
a2b2 = beta^2-alpha^2;
for i = 2:N-1
    ab2 = 2*i+alpha+beta;
    a(i) = a2b2/((ab2-2)*ab2);
    b(i) = sqrt((4*i*(i+alpha)*(i+beta)*(i+alpha+beta))/((ab2^2-1)*ab2^2));
end
ab2 = 2*N+alpha+beta;
a(N) = a2b2/((ab2-2)*ab2);

% Construct tridiagonal matrix
Matrix = diag(a)+diag(b,1)+diag(b,-1);

% Use eigenvalues and eigenvectors of matrix to compute nodes and weights
[evec,eval] = eig(Matrix);
weights = evec(1,:);
[nodes,index] = sort(diag(eval));
weights = weights(index);
nodes = nodes';
weights = 2^(alpha+beta+1)*gamfun(alpha+1)*gamfun(beta+1)/ ...
    gamfun(alpha+beta+2).*weights.^2;
