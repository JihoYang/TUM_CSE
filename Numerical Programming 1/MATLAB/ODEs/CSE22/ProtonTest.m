close all;

[t,x] = ode15s(@proton,[0 8e5],[0;1;0]);

semilogx(t,x(:,3),'b.','MarkerSize',7); shg;
xlim([1e-16,1e8]); ylim([0,3.5e-17]);

pause; clf;

options = odeset('AbsTol',1e-20);

[t,x] = ode15s(@proton,[0 8e5],[0;1;0],options);

semilogx(t,x(:,3),'b.','MarkerSize',7); shg;
xlim([1e-16,1e8]); ylim([0,3.5e-17]);

pause; clf;

options = odeset('AbsTol',1e-20,'InitialStep',1e-20);

[t,x] = ode15s(@proton,[0 8e5],[0;1;0],options);

semilogx(t,x(:,3),'b.','MarkerSize',7); shg;
xlim([1e-16,1e8]); ylim([0,3.5e-17]);
