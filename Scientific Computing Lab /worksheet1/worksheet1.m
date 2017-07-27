%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Worksheet 1       %%
%% Scientific Computing Lab %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------%
%        Written by:       %
%--------------------------%
%        Nathan Brei       %
%         Jiho Yang        %
%   Tord Kriznik Sørensen %
%--------------------------%


close all;
clear all;
clc;


%--------------------%
%-- Initialisation --%
%--------------------%


%pd:   Derivative of p.
%p0:   Initial value of p0 at t0.
%t0:   Initial time.
%dt:   Timestep.
%tmax: End time.
%tt:   Time to be used in exact solution.
%pp:   Exact solution to the ODE given by pd.

% We included t in the function pd because we wanted a general method which can
% take in functions explicitly dependant on t for later use/testing.

pd   = @(t,p) (1-p/10)*p;
p0   = 1;
t0   = 0;
dt = [1/8 1/4 1/2 1];
tmax = 5;
tt = [];
pp = []; 

% E: Matrix to analyse the errors of the different methods.
% E(1,:) Contains different values of dt
% E(2,:) Contains the approximation error using p_exact
% E(3,:) Contains the error augmentation when dt = 2*dt. (First entry not applicable)
% E(4,:) Contains the approximation error p_bestapprox

E_euler = zeros(4);
E_euler(1,:) = dt;

E_heun = zeros(4);
E_heun(1,:) = dt;

E_rk = zeros(4);
E_rk(1,:) = dt;


% p: The approximated values by the different methods for different dt
% t: The timesteps for the different methods for different dt

p_euler = cell(1,4);
p_heun  = cell(1,4);
p_rk    = cell(1,4);

t_euler = cell(1,3);
t_heun  = cell(1,3);
t_rk    = cell(1,3);


%--------------------------%
%-- Solution Computation --%
%--------------------------%


for i=1:length(dt)

    % No timestep vector (t) is defined in the main script, but instead it
    % is calculated inside the functions and returned.
    % This is not necessary here since all the t vectors used for a given
    % dt are the same for all methods. However, we wanted each function to
    % be self contained, and hence for the completeness of each function,
    % the t vectors are computed within them.
    
    % We return the timesteps (t) because we want to use the same
    % convention as ode45
    
    % Calculating the approximated values for p.
    [t_euler{i}, p_euler{i}] = explicit_euler(pd,p0,t0,dt(i),tmax);
    [t_heun{i}, p_heun{i}]   = heun(pd,p0,t0,dt(i),tmax);
    [t_rk{i},p_rk{i}]        = rungekutta(pd,p0,t0,dt(i),tmax);

    % Calculating the exact values.
    tt{i} = [t0:dt(i):tmax];
    pp{i} = 10./(1+9*exp(-tt{i}));

    % Finding the real error.
    E_euler(2,i) = getError(p_euler{i},pp{i},dt(i));
    E_heun(2,i)  = getError(p_heun{i},pp{i},dt(i));
    E_rk(2,i)    = getError(p_rk{i},pp{i},dt(i));

    % Comparison of errors.
    if  i>1
        E_euler(3,i) = E_euler(2,i)/ E_euler(2,i-1);
        E_heun(3,i)  = E_heun(2,i)/ E_heun(2,i-1);
        E_rk(3,i)    = E_rk(2,i)/ E_rk(2,i-1);
    end

    % Finding error based on best approximation.
    % Prune the 'best' results y_*{1} by taking every nth value
    % so we can compare directly against y_*{i}

    scale = dt(i)/dt(1);
    E_euler(4,i) = getError(p_euler{i}, p_euler{1}(1:scale:end), dt(i));
    E_heun(4,i)  = getError(p_heun{i},  p_heun{1}(1:scale:end),  dt(i));
    E_rk(4,i)    = getError(p_rk{i},    p_rk{1}(1:scale:end),    dt(i));

    % Sanity-check: Approximate error should generally be less than
    % absolute error, because the 'best' approximate result should lie
    % between the exact result and any other approximate result, as long
    % as we encounter no instabilities.
    
end


%------------%
%-- Output --%
%------------%


fprintf('<strong>Explanation of tables </strong>')
disp(' ')
disp(' ')
disp('First row:  Timestep.')
disp('Second row: The approximation error using the exact solution for p as reference.')
disp('Third row:  The increase in error from doubling the timestep.')
disp('Fourth row: The approximation error using the best approximated solution as reference.')

disp(' ')
disp(' ')
fprintf('<strong>Results </strong>')
disp(' ')

E_euler
E_heun
E_rk   


%-----------%
%-- Plots --%
%-----------%

figure (1);
title('Explicit Euler Method for Different Timesteps (dt)');
hold on;

plot(tt{1},pp{1},'b','LineWidth',1.5);
plot(t_euler{1}, p_euler{1},'ro--','LineWidth',1.5);
plot(t_euler{2}, p_euler{2},'gx--','LineWidth',1.5);
plot(t_euler{3}, p_euler{3},'y+--','LineWidth',1.5);
plot(t_euler{4}, p_euler{4},'c*--','LineWidth',1.5);

set(gca, 'FontSize', 15)
xlabel('t','FontSize',20);
ylabel('P(t)','FontSize',20);
legend('Analytical Solution','Explicit Euler (dt=1/8)', 'Explicit Euler (dt=1/4)','Explicit Euler (dt=1/2)', 'Explicit Euler (dt=1)','Location','northwest');

box on; grid on;

figure (2);
title('Heun''s Method for Different Timesteps (dt)');
hold on;

plot(tt{1},pp{1},'b','LineWidth',1.5);
plot(t_heun{1}, p_heun{1},'ro--','LineWidth',1.5);
plot(t_heun{2}, p_heun{2},'gx--','LineWidth',1.5);
plot(t_heun{3}, p_heun{3},'y+--','LineWidth',1.5);
plot(t_heun{4}, p_heun{4},'c*--','LineWidth',1.5);

set(gca, 'FontSize', 15)
xlabel('t','FontSize',20);
ylabel('P(t)','FontSize',20);
legend('Analytical Solution','Heun (dt=1/8)', 'Heun (dt=1/4)','Heun (dt=1/2)', 'Heun (dt=1)', 'Location','northwest');

box on; grid on;

figure (3);
title('Runge Kutta Method for Different Timesteps (dt)');
hold on;

plot(tt{1},pp{1},'b','LineWidth',1.5);
plot(t_rk{1}, p_rk{1},'ro--','LineWidth',1.5);
plot(t_rk{2}, p_rk{2},'gx--','LineWidth',1.5);
plot(t_rk{3}, p_rk{3},'y+--','LineWidth',1.5);
plot(t_rk{4}, p_rk{4},'c*--','LineWidth',1.5);

set(gca, 'FontSize', 15)
xlabel('t','FontSize',20);
ylabel('P(t)','FontSize',20);
legend('Analytical Solution','Runge Kutta (dt=1/8)', 'Runge Kutta (dt=1/4)','Runge Kutta (dt=1/2)', 'Runge Kutta (dt=1)', 'Location','northwest');

box on; grid on;


