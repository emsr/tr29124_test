clear all;close all;
tiny=0.06447*eps;
y=logspace(-20.,4,71);
x=linspace(-200,200,40001);
[x,y]=meshgrid(x,y);
z=complex(x,y);

%  tic, operation, toc: prints the number of seconds required for the
%  operation

tic; [w1]=Faddeyeva(z,tiny);toc  % Present algorithm
% tic; [w1]=Faddeyeva(z);toc       % Present algorithm
% tic; [w2]=cerf(z);toc            % Hui et al [1978]
% tic; [w3]=humlicek_rev(z);toc    % Humlicek  [1982], original
% tic; [w3]=w4(z);toc              % Humlicek  [1982], modified
% tic; [w4]=w(z);toc               % Poppe & Wijers [1990]
% tic; [w5]=Weideman(z,128);toc    % Weidemann  [1994]

% Calculations of partial derivatives
dVdx=-2*real(z.*w1);             % Partial derivative of V w.r.t. x
dVdy=2*imag(z.*w1)-2/sqrt(pi);   % Partial derivative of V w.r.t. y