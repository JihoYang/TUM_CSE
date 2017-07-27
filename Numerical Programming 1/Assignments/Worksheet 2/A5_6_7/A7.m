function A7(n)
x = cos(pi/n*(0:n));
y = 1./(1+25*x.^2);
lambda = (-1).^(0:n); lambda([1,end]) = lambda([1,end])/2;
xx = linspace(-1,1,1e4); yy = 1./(1+25*xx.^2);
num = 0; den = 0;
for j=1:n+1
    num = num + lambda(j)*y(j)./(xx-x(j));
    den = den + lambda(j)./(xx-x(j));
end
p = num./den; 

P = polyfit(x,y,n);


%plot(xx,yy,xx,p,x,y,'*',xx,polyval(P,xx)); 
plot(xx,yy,xx,p,x,y,'*');

legend({'function','barycentric formula','interp. points','poly-fit/val'});
shg;
end