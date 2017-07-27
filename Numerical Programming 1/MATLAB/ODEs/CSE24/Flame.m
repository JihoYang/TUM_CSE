f = @(t,y) y^2*(1-y);

epsilon = 1e-5;

disp(sprintf('ode45:\n'));

disp(sprintf('prediced number of steps = %i',ceil(1/epsilon/3.31)));

options = odeset('Refine',1,'Stats','on');
tic, [t,y] = ode45(f,[0,2/epsilon],epsilon,options); toc

figure(1)
plot(t,y,'bo');

disp(sprintf('\n ode23s:\n'));
options = odeset('Refine',1,'Stats','on');
tic, [t,y] = ode23s(f,[0,2/epsilon],epsilon,options); toc

figure(2)
plot(t,y,'bo');