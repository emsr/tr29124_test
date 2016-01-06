function [f] = gamfun_godfrey(z)
% function [f] = gamfun_godfrey(z)
% 
% Computes the Gamma function using Godfrey's method, involving a Lanczos
% series.
% 
% Based on work of Paul Godfrey, discussed by Gerard Michon at:
%    http://www.numericana.com/answer/info/godfrey.htm
% 
% Copyright John W. Pearson 2014


% Parameter choice of Godfrey
g = 607/128;

% Vector of coefficients of the Lanczos series
c = [0.99999999999999709182; 57.156235665862923517; -59.597960355475491248; ...
    14.136097974741747174; -0.49191381609762019978; 3.3994649984811888699e-5; ...
    4.6523628927048575665e-5; -9.8374475304879564677e-5; 1.5808870322491248884e-4; ...
    -2.1026444172410488319e-4; 2.1743961811521264320e-4; -1.6431810653676389022e-4; ...
    8.4418223983852743293e-5; -2.6190838401581408670e-5; 3.6899182659531622704e-6];

% Compute parameter zp
zp = (z+g-0.5)^(0.5*(z-0.5));

% Sum series
s = 0;
for j = length(c)-1:-1:1
    s = s+c(j+1)/(z+j-1);
end

% Use series sum to compute Gamma function
f = (sqrt(2*pi)*(c(1)+s))*((zp*exp(-(z+g-0.5)))*zp);

% Case where z is real
if imag(z)==0
    if z==fix(z) && real(z)<=0 % z is a negative integer
        f = Inf;
    else
        f = real(f); % if z is real, then Gamma function is real
    end
end
