function [y,err,nfcn] = quad24(f,a,b,tol)
t = a; y = 0; 
h = (b-a)/100;                              % initial step size
f_left = f(a); nfcn = 1; err=0;
while t < b
    f_middle = f(t+h/2); f_right = f(t+h);  % evaluate f
    nfcn=nfcn+2;                            % count f evaluations   
    Q1 = h*(f_left+f_right)/2;              % trapezoidal rule (order 2)
    Q2 = h*(f_left+4*f_middle+f_right)/6;   % Simpson rule     (order 4)
    delta = abs(Q1-Q2);                     % error estimate 
    if delta < tol
        t = t + h;
        y = y + Q2; err=err+delta;
        f_left = f_right;
    end
    h = min(0.9*h*(tol/delta)^(1/3),b-t);   % step size control
end 
