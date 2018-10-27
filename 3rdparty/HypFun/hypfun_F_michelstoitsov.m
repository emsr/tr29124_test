function [h] = hypfun_F_michelstoitsov(a,b,c,z,tol)
% function [h] = hypfun_F_michelstoitsov(a,b,c,z,tol)
% 
% Computes the hypergeometric function \mathbf{F}(a,b;c;z), for z near the
% unit disc, using Michel and Stoitsov's method.
% 
% Copyright John W. Pearson 2014


% Compute coefficient as recommended in paper of Michel and Stoitsov
if abs(z)<=1
    r0 = 0.9;
elseif abs(z)>1
    r0 = 1.1;
end

% Compute quantity z_0
z0 = r0*z/abs(z);

% Start off sum
qmid = gamfun(c)*hypergeom([a,b],c,z0);
qnew = a*b/c*gamfun(c)*hypergeom([a+1,b+1],c+1,z0);
a1 = z-z0;
h = qmid+qnew*a1;

% Compute sum
for j = 2:500
    a1 = a1*(z-z0);
    qold = qmid; qmid = qnew;
    qnew = 1/(z0*(1-z0)*j)*(((j-2)*(2*z0-1)-c+(a+b+1)*z0)*qmid+ ...
        (a+j-2)*(b+j-2)/(j-1)*qold);
    h = h+qnew*a1;
    if abs(qnew*a1)/abs(h)<tol
        break
    end
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% Compute \mathbf{F}
h = h/gamfun(c);
