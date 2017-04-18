function w = Faddeyeva(z,tiny)
%----------
%   w = Faddeyeva(z,tiny) is the Faddeyeva or the plasma dispersion
%   function, w(z) = exp(-z^2) * erfc(-i*z) where z=x+iy and erfc(z) is the
%   complex  complementary error function of z; z is usually an array (with
%   one or two dimensions) but can be a single scalar as well.
%----------
%   The parameter "tiny" can be assigned values between the default value
%   tiny=1e-4 (for lowest accuracy & shortest run time) and tiny_min<eps
%   (for highest accuracy with longer run time) where eps is the floating
%   point relative accuracy in the computational platform. However, it is
%   not recommended for tiny_min to be << eps. A default value of tiny_min=
%   0.06447*eps is used herein and any smaller input value for
%   tiny will be automatically changed to this value with a warning
%   message.
%-------------
%   The accompanying driver code Faddeyeva_driver.m can be run for 
%   computation of the Faddeyeva function and the partial derivatives
%   of its real part, V(x,y) on an array of the complex variable z.The 
%   partial derivatives of the imaginary part, L(x,y) are given simply 
%   by Eq. (23) in the manuscript. 
%   An example of generating an array of z is included in the driver code.  
%-------------
%   Authors: M. R. Zaghloul & A. N. Ali,
%   United Arab Emirates University, May 28,2011
%----------

% Machine-related parameters
Rmin=realmin;                     % Smallest positive floating point number
eps0=eps;                         % Floating point relative accuracy
sqrt_log_Rmin=sqrt(-log(Rmin));

%-----------
% Check inputs, at least one input argument is needed
if nargin<1
    error('Faddeyeva.m needs at least one input argument');
    w=NaN+i*NaN;
    return
elseif nargin==1, tiny=[];
    % Set the default value of the parameter "tiny"
    if isempty(tiny),  tiny=1e-4; end
elseif nargin==2
    % Set the maximum & minimum values of the parameter tiny with warning
    % messages when the input tiny goes beyond the limits
    if ~isempty(tiny) && tiny > 1e-4;
        tiny=1e-4;
        ws1=sprintf(['tiny must be <= 1e-4.'...
            'The default value of tiny=1e-4 has been used']);
        warning(ws1);
    elseif nargin==2 && ~isempty(tiny) && tiny < (0.06447*eps0);
        tiny = (0.06447*eps0);
        ws2=sprintf(['tiny is less than tiny_miniumum=0.06447*eps \n'...
            'The value of tiny_minimum has been used         ']);
        warning(ws2);
    end
end

% For purely imaginary input values, use the built-in function erfcx.
if isreal(-1i*z)
    w=erfcx(imag(z));
    return
end

% Faddeyeva cannot calculate w for y negative & exp(y^2-x^2)>=the largest
% floating point number due to unavoidable over flow problems
if ~isempty(find(imag(z)<0, 1));
    neg_y_idx=find((imag(z).^2-real(z).^2)>=sqrt_log_Rmin^2, 1);
    if ~isempty(neg_y_idx)
        es=sprintf(['Input array, z, has points with imag(z)<0 & \n',...
            'exp(imag(z)^2-real(z)^2)>=the largest floating point number.\n', ...
            'Faddeyeva cannot calculate w for these points ',...
            'due to unavoidable overflow problems ']);
        error(es)
    end
end
%------------
% Calculate the value of the parameter "a" corresponding to "tiny"

a=realsqrt(-(pi*pi)/reallog(tiny/2));          % parameter a in Eq. (8)
tiny=max(tiny,eps0);

%-------------------

% Repeatedly used parameters & constants
half_a=.5*a;
a_sqr=a*a;
two_a=2*a;
two_a_sqr=2*a_sqr;
four_a_sqr=2*two_a_sqr;
a_pi=a/pi;
two_a_pi=2*a_pi;
one_sqrt_pi=(1/sqrt(pi));

% Reset outputs
sizein=size(z);
w=repmat(NaN,sizein);
idxmax=prod(sizein);
%--------------------

% Calculate erfcx(y).
erfcx_y=erfcx(abs(imag(z)));
%--------------------
for idx=1:idxmax;
    x=real(z(idx));
    y=imag(z(idx));
    xsign=sign(x);
    ysign=sign(y);
    x=xsign*x;
    y=max(Rmin,ysign*y);
    
    x_sqr=x*x;
    y_sqr=y*y;
    two_yx=2*y*x;
    two_a_x=two_a*x;
    exp_x_sqr=exp(-x_sqr);
    cos_2yx=cos(two_yx);
    sin_2yx=sin(two_yx);
    
    % For x=0, use the asymptotic expressions for x--->0 from Eq. (6) in
    % the manuscript
    if x==0
        w(idx)=erfcx_y(idx)*(1+1i*xsign*(x/y));
    else
        
        V_old=exp_x_sqr*(erfcx_y(idx)*cos_2yx+(two_a_pi/y)*sin(two_yx/2)^2);
        L_old=(-erfcx_y(idx)+a_pi/y);
        
        % Initialization of the sums Sigma3, Sigma5 & Sigma4_5
        Sigma3=Rmin;
        Sigma5=Rmin;
        Sigma4_5=0;
        
        delta3=1;
        delta5=1;
        
        n=0;
        n3=ceil(x/a);
        n3_3=n3-1;
        
        if (sqrt_log_Rmin-x)>0
            
            % Initialization of the sums Sigma1, Sigma2 & Sigma4
            Sigma1=Rmin;
            Sigma2=Rmin;
            Sigma4=Rmin;
            
            exp1=exp(-two_a_x);
            exp2=exp((four_a_sqr*n3-2*two_a_x)-two_a_sqr);
            exp3=exp(-((two_a_sqr*n3-two_a_x)-two_a_sqr));
            del2_tmp=1.0;
            del3_tmp=exp(-((a_sqr*n3*n3-two_a_x*n3-two_a_sqr*n3)+x_sqr+two_a_x+a_sqr));
            del3_3_tmp=exp(a_sqr-(two_a_sqr*n3-two_a_x));
            
            while (delta3>=tiny || delta5>=tiny && n<=50)
                n=n+1;
                den1=a_sqr*n*n+y_sqr;
                exp_del1=exp(-(a_sqr*n*n))/den1;
                del2_tmp=del2_tmp*exp1;
                
                if n3_3>=1
                    del3_tmp=del3_tmp*exp3;
                    exp3_den=del3_tmp*exp_del1*(den1/(a_sqr*n3*n3+y_sqr));
                    del3_3_tmp=del3_3_tmp*exp2;
                    exp3_3_den=exp3_den*del3_3_tmp*((a_sqr*n3*n3+y_sqr)/(a_sqr*n3_3*n3_3+y_sqr));
                    del5=(n3_3*exp3_3_den+n3*exp3_den);
                    del3=(exp3_3_den+exp3_den);
                else
                    del3_tmp=del3_tmp*exp3;
                    del3=del3_tmp*exp_del1*(den1/(a_sqr*n3*n3+y_sqr));
                    del5=(n3*del3);
                end
                
                delta3=del3/Sigma3;
                delta5=del5/Sigma5;
                
                Sigma1=Sigma1+exp_del1;
                Sigma2=Sigma2+del2_tmp*exp_x_sqr*exp_del1;
                Sigma3=Sigma3+del3;
                Sigma4=Sigma4+n*del2_tmp*exp_x_sqr*exp_del1;
                Sigma5=Sigma5+del5 ;
                
                if x>=5e-4
                    Sigma4_5=-Sigma4+Sigma5;
                else
                    Sigma4_5=Sigma4_5+2*n*n*two_a_x*exp_x_sqr*exp_del1*(1+1.666666666666667e-001*(two_a_x*n)^2+...
                        8.333333333333333e-003*(two_a_x*n)^4);
                end
                
                n3=n3+1;
                n3_3=n3_3-1;
                
            end
            
            if y <= 5e0 && two_yx>Rmin
                w(idx)=V_old+y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+.5*(Sigma2+Sigma3))+...
                    1i* xsign*(sin_2yx*exp_x_sqr*(L_old+two_a_pi*y*Sigma1)+two_a_pi*half_a*Sigma4_5);
            elseif y <= 5e0 && two_yx<=Rmin
                w(idx)=V_old+y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+.5*(Sigma2+Sigma3))+...
                    1i* xsign*(2*y*exp_x_sqr*(x*L_old+x*two_a_pi*y*Sigma1)+two_a_pi*half_a*Sigma4_5);
            else
                w(idx)=V_old+y*two_a_pi*(-cos_2yx*exp_x_sqr*Sigma1+.5*(Sigma2+Sigma3))+...
                    1i* xsign*(sin_2yx*exp_x_sqr*min(0,abs(L_old+(two_a_pi*y*Sigma1)))+two_a_pi*half_a*Sigma4_5);
            end
            
        elseif x>=sqrt_log_Rmin && x<1e15
            exp2=exp((four_a_sqr*n3-2*two_a_x)-two_a_sqr);
            del3_3_tmp=exp(a_sqr+(two_a_x-two_a_sqr*n3));
            
            while (delta3>=tiny || delta5>=tiny && n<=50)
                n=n+1;
                
                if n3_3>=1
                    exp3_den=exp(-(a*n3-x)*(a*n3-x))/(a_sqr*n3*n3+y_sqr);
                    del3_3_tmp=del3_3_tmp*exp2;
                    exp3_3_den=exp3_den*del3_3_tmp*((a_sqr*n3*n3+y_sqr)/(a_sqr*n3_3*n3_3+y_sqr));
                    del5=(n3_3*exp3_3_den+n3*exp3_den);
                    del3=(exp3_3_den+exp3_den);
                else
                    del3=exp(-(a*n3-x)^2)/(a_sqr*n3*n3+y_sqr);
                    del5=n3*del3;
                end
                
                delta3=del3/Sigma3;
                delta5=del5/Sigma5;
                Sigma3=Sigma3+del3 ;
                Sigma5=Sigma5+del5;
                n3=n3+1;
                n3_3=n3_3-1;
                
            end
            w(idx)=V_old+y*a_pi*Sigma3+1i*xsign*(sin_2yx*exp_x_sqr*L_old+two_a_pi*half_a*Sigma5);
        else
            w(idx)=one_sqrt_pi*((y+1i*xsign*x)/(x_sqr+y_sqr)) ;
        end
    end
    
    if ysign<0
        two_exp_x_sqr_ysqr=2*exp(-x_sqr+y_sqr);
        w(idx)=two_exp_x_sqr_ysqr*cos_2yx-real(w(idx))-1i*(-xsign*two_exp_x_sqr_ysqr*sin_2yx-imag(w(idx)));
    end
end

