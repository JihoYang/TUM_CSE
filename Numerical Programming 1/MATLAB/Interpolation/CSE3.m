x = linspace(-1,1,7)

x =

  Columns 1 through 2

                        -1        -0.666666666666667

  Columns 3 through 4

        -0.333333333333333                         0

  Columns 5 through 6

         0.333333333333333         0.666666666666667

  Column 7

                         1

format short g
x

x =

  Columns 1 through 5

           -1     -0.66667     -0.33333            0      0.33333

  Columns 6 through 7

      0.66667            1

y = [0 0 0 0 0 1 0];
p = polyfit(x,y,6);
xx = linspace(-1,1,100); yy = polyval(p,xx);
plot(xx,yy,x,y,'r*')

x = linspace(-1,1,21); 
y = zeros(size(x)); y(19)=1;
p = polyfit(x,y,20);
xx = linspace(-1,1,250); yy = polyval(p,xx);
plot(xx,yy,x,y,’r*’); shg

n = 1000;
2/pi*log(n+1)+1

ans =

       5.3982

n = 8; x = cos(pi/n*(0:n)); plot(x,zeros(size(x)),’r*’); shg
n = 50; x = cos(pi/n*(0:n)); plot(x,zeros(size(x)),’r*’); shg

n = 8;
x = cos(pi/n*(0:n));
y = 1./(1+25*x.^2);
for j=1:n+1, lambda(j) = prod(1./(x(j)-x([1:j-1,j+1:end]))); end
xx = linspace(-1,1,100); yy = 1./(1+25*xx.^2);
omega = 1; s = 0;
for j=1:n+1, omega = omega.*(xx-x(j)); s = s + y(j)*lambda(j)./(xx-x(j)); end
p = omega.*s; plot(xx,yy,xx,p,x,y,'*')

ChebyshevInterpolation
ChebyshevInterpolation2
