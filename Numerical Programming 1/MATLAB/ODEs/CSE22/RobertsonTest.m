close all;

h = 1e-5;
options = odeset('InitialStep',h);

y0 = [1;0;0]; 

T = 4e6; 
tic, [t,y] = ode15s(@Robertson,[0,T],y0,options); toc;

y(:,2) = 1e4*y(:,2);

semilogx(t,y,'.','MarkerSize',7); shg;
xlim([h,T]);

pause; clf;

T = 4e10; 
tic, [t,y] = ode15s(@Robertson,[0,T],y0,options); toc;

y(:,2) = 1e4*y(:,2);

semilogx(t,y,'.','MarkerSize',7); shg;
xlim([h,T]);

pause; clf;

T = 4e13; 
tic, [t,y] = ode15s(@Robertson,[0,T],y0,options); toc;

y(:,2) = 1e4*y(:,2);

semilogx(t,y,'.','MarkerSize',7); shg;
xlim([h,T]);

pause; clf;

T = 4e13; 
options = odeset(options,'NonNegative',[1 2 3]);
tic, [t,y] = ode15s(@Robertson,[0,T],y0,options); toc;

y(:,2) = 1e4*y(:,2);

semilogx(t,y,'.','MarkerSize',7); shg;
xlim([h,T]);