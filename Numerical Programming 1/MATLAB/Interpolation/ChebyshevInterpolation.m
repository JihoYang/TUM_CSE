% first barycentric formula with Chebyshev nodes
n = 50;
x = cos(pi/n*(0:n));
% x = linspace(-1,1,n+1); % if you want to see the equistant desaster
y = 1./(1+25*x.^2);
for j=1:n+1, lambda(j) = prod(1./(x(j)-x([1:j-1,j+1:end]))); end
xx = linspace(-1,1,100); yy = 1./(1+25*xx.^2);
omega = 1; s = 0;
for j=1:n+1, omega = omega.*(xx-x(j)); s = s + y(j)*lambda(j)./(xx-x(j)); end
p = omega.*s; plot(xx,yy,xx,p,x,y,'*'); shg;