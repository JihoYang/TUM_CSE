function res = F_func(xx)

res=xx;
x = xx(1); y = xx(2);
x2 = x^2;
y2 = y^2;
res(1) = (1 - 4*x2)*(x + 2*y - 2*y^3) + 4*x*y2*(x2 - 1) + 3*x*y^4;
res(2) = 2*x*(1 - x2)^2 + y*(1 - 2*x2 + 4*x^4) + 6*x*y2*(x2 - y2) - 3*y^5;

end
