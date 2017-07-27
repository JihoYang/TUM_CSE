function dx=threebody(t,x)

dx = x;

mu = 0.012277471; ms = 1 - mu;
r1 = norm(x(1:2)-[-mu;0]);
r2 = norm(x(1:2)-[ ms;0]);

dx(1:2) = x(3:4);
dx(3) = x(1) + 2*x(4) - ms*(x(1)+mu)/r1^3 - mu*(x(1)-ms)/r2^3;
dx(4) = x(2) - 2*x(3) - ms*     x(2)/r1^3 - mu*     x(2)/r2^3;