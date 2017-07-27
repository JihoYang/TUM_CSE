function [y,err,nfcn] = A14(f,a,b,tol,Iscal)
% adaptive Simpson quadrature with embedded trapezoidal rule
% function [y,err,nfcn] = A14(f,a,b,tol,Iscal)
%
% Remark: This function uses only the relative error concept.
%         There is no AbsTol (i.e. AbsTol=0)!!
if nargin<5, Iscal=1; end
t = a; y = 0; tolIscal = tol*Iscal; 
h = (b-a)*min(1/100,tol);                   % initial step size
f_left = f(a); nfcn=1; err=0;
while t < b
    f_middle = f(t+h/2); f_right = f(t+h);  % evaluate f
    nfcn=nfcn+2;                            % count f evaluations
    Q1 = h*(f_left+f_right)/2;              % trapezoidal rule (order 2)
    Q2 = h*(f_left+4*f_middle+f_right)/6;   % Simpson rule     (order 4)
    delta = abs(Q1-Q2);                     % error estimate for trapez
    epsilon = (delta^2)/max(abs(Q1),Iscal*eps);  % error estimate for Simpson
    if epsilon < tolIscal
        t = t + h;
        f_left = f_right;                   % reuse f value next time
        y = y + Q2; err=err+epsilon;       
    end
    h = min(min(0.9*h*(tolIscal/epsilon)^(1/5),b-t),3*h); % step size control
end 
