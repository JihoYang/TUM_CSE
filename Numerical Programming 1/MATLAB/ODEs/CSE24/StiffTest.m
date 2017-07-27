function StiffTest(tol)

close all;

%tol = 1e-2;

options = odeset('Refine',1,'RelTol',tol,'AbsTol',1e-6*tol,'Stats','on');

y0 = [1;0;0]; T=0.3;

tic, [t,y] = ode45(@Robertson,[0,T],y0,options); toc;

y(:,2) = 1e4*y(:,2);

plot(t,y,'.','MarkerSize',7); shg;
xlim([0,T]);
end