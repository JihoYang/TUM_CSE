function [y,n] = quad24(f,a,b,tol)

t = a; y = 0; n = 1;
h = (b-a)/100;                              % initial step size
while t < b
    Q1 = h*(f(t)+f(t+h))/2;                 % trapezoidal rule (order 2)
    Q2 = h*(f(t)+4*f(t+h/2)+f(t+h))/6;      % Simpson rule     (order 4)
    delta = abs(Q1-Q2);                     % error estimate 
    if delta < tol
        t = t + h;
        y = y + Q2;
    end
    h = min(0.9*h*(tol/delta)^(1/3),b-t);   % step size control
    n = n + 1;
end 
