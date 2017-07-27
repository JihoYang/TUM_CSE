f = @(x) cosh(sin(x));

% non-periodic case

val = 2.283194520910577988; % integral from x=-1 to x=1
nmax = 500;
for n=2:nmax
    x = linspace(-1,1,n);   
    int = trapezoidal(f,x);
    err(n) = abs(int - val);
end

subplot(1,2,1)
loglog(2:nmax,err(2:nmax),'b.','MarkerSize',10)
grid on

% periodic case

val = 7.954926521012845275; % integral from x=-pi to x=pi
nmax = 20;
for n=2:nmax
    x = linspace(-pi,pi,n);   
    int = trapezoidal(f,x);
    err(n) = abs(int - val);
end

subplot(1,2,2)
semilogy(2:nmax,err(2:nmax),'b.','MarkerSize',15)
grid on
    