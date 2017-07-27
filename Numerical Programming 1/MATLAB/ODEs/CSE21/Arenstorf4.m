close all;

T = 17.06521656015796; mu = 0.0123;
y0 = [0.994; 0; 0; -2.001585106379082];

tic, [t,x] = Euler(@threebody,linspace(0,T,100000),y0); toc;

hold on;
plot(x(:,1),x(:,2),'b.','MarkerSize',7); 
plot(-mu,0,'r.','MarkerSize',10); plot(1-mu,0,'r.','MarkerSize',10); shg;
hold off;

pause; clf;

options = odeset('RelTol',1e-5,'Refine',1);
tic, [t,x] = ode45(@threebody,[0,T],y0,options); toc

hold on;
plot(x(:,1),x(:,2),'b.','MarkerSize',7); 
plot(-mu,0,'r.','MarkerSize',10); plot(1-mu,0,'r.','MarkerSize',10); shg;
hold off;

pause; clf;

tic, sol = ode45(@threebody,[0,T],y0,options); toc

tt = linspace(0,T,10000);
xx = deval(sol,tt)';

hold on;
plot(xx(:,1),xx(:,2),'g-','LineWidth',2);
plot(x(:,1),x(:,2),'b.','MarkerSize',7); 
plot(-mu,0,'r.','MarkerSize',20); plot(1-mu,0,'r.','MarkerSize',20); shg;
hold off;