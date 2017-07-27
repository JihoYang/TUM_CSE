close all;

T = 6.192169331; mu = 0.0123;
y0 = [1.2; 0; 0; -1.049357510];

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
plot(xx(:,1),xx(:,2),'g-','LineWidth',2,'LineSmoothing','on');
plot(x(:,1),x(:,2),'b.','MarkerSize',7); 
plot(-mu,0,'r.','MarkerSize',20); plot(1-mu,0,'r.','MarkerSize',20); shg;
hold off;