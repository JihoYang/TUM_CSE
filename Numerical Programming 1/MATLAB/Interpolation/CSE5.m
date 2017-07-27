n = 10;
f = @(x) x.^2 + 10./(sin(x)+1.2);
x = linspace(-2,8,n); y = f(x);
xx = linspace(-2,8,250); yy = f(xx);
p = polyfit(x,y,n-1); pp = polyval(p,xx);
plot(xx,yy,'b-',xx,pp,'r-',x,y,'r*'), shg

s = interp1(x,y,xx); 
plot(xx,yy,'b-',xx,s,'r-',x,y,'r*'), shg

s = interp1(x,y,xx,'spline');
figure(2)
plot(xx,yy,'b-',xx,s,'r-',x,y,'r*'), shg

plot(xx,yy,'b-',xx,s,'r-'), shg

n = 100;
x = linspace(-2,8,n); y = f(x);
s = interp1(x,y,xx,'spline');
tic, s = interp1(x,y,xx,'spline'); toc
Elapsed time is 0.001062 seconds.

plot(xx,yy,'b-',xx,s,'r-'), shg
p = polyfit(x,y,n-1); pp = polyval(p,xx);
[Warning: Polynomial is badly conditioned. Add points with distinct X
values, reduce the degree of the polynomial, or try centering and
scaling as described in HELP POLYFIT.] 

plot(xx,yy,'b-',xx,pp,'r-',x,y,'r*'), shg
plot(xx,yy,'b-',xx,s,'r-'), shg

n = 10000;
x = linspace(-2,8,n); y = f(x);
xx = linspace(-2,8,250); yy = f(xx);
tic, s = interp1(x,y,xx,'spline'); toc
Elapsed time is 0.003183 seconds.

n = 1e6;
x = linspace(-2,8,n); y = f(x);
tic, s = interp1(x,y,xx,'spline'); toc
Elapsed time is 0.140373 seconds.
