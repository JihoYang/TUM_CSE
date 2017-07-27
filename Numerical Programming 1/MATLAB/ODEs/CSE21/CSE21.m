% Arenstorf4
% Arenstorf3
% Arenstorf2

f = @(t,y) 5*(y-t^2);
y_exact = @(t) t^2 + 0.4*t + 0.08;
T = 10;
[t,y] = ode45(f,[0,T],0.08);
[y(end) y_exact(T)]

options = odeset('RelTol',5e-14);
[t,y] = ode45(f,[0,T],0.08,options);
[y(end) y_exact(T)]

[t,y] = ode45(f,[20,T],1e8,options);
[y(end) y_exact(T)]