%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Worksheet 2       %%
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
%err:  Accuracy limit for the Newton iteration in each time step.

% We included t in the function pd because we wanted a general method which can
% take in functions explicitly dependant on t for later use/testing.

pd   = @(t,p) 7.*(1-p/10).*p;
pdd  = @(t,p) 7.*(1-p/5);
p0   = 20;
t0   = 0;
dt = [1/32 1/16 1/8 1/4 1/2];
tmax = 5;
tt = [];
pp = []; 

% E: Matrix to analyse the errors of the different methods.
% E(1,:) Contains different values of dt
% E(2,:) Contains the approximation error using p_exact
% E(3,:) Contains the error augmentation when dt = 2*dt. (First entry not applicable)
% E(4,:) Contains the approximation error p_bestapprox

rows = 3; %4 %when calculating approximate errors
E_euler = zeros(rows,length(dt));
E_euler(1,:) = dt;

E_heun = zeros(rows, length(dt));
E_heun(1,:) = dt;

E_implicit_euler = zeros(rows, length(dt));
E_implicit_euler(1,:) = dt;

E_adams_moulton = zeros(rows, length(dt));
E_adams_moulton(1,:) = dt;

E_adams_moulton_linear1 = zeros(rows, length(dt));
E_adams_moulton_linear1(1,:) = dt;

E_adams_moulton_linear2 = zeros(rows, length(dt));
E_adams_moulton_linear2(1,:) = dt;

% p: The approximated values by the different methods for different dt
% t: The timesteps for the different methods for different dt

p_euler = cell(1,length(dt));
p_heun  = cell(1,length(dt));
p_implicit_euler = cell(1,length(dt));
p_adams_moulton = cell(1,length(dt));
p_adams_moulton_linear1 = cell(1,length(dt));
p_adams_moulton_linear2 = cell(1,length(dt));

t_euler = cell(1,length(dt));
t_heun  = cell(1,length(dt));
t_implicit_euler = cell(1,length(dt));
t_adams_moulton = cell(1,length(dt));
t_adams_moulton_linear1 = cell(1,length(dt));
t_adams_moulton_linear2 = cell(1,length(dt));

%TODO: Print stable in a good way.
% stable: Cell array of boolean values. which show wheter each
%         function is stable for given dt.
%stable = cell(length(dt) + 1,1);
%stable{1} = dt;

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
    [t_implicit_euler{i}, p_implicit_euler{i}] = implicit_euler(pd,pdd,p0,t0,dt(i),tmax);
    [t_adams_moulton{i}, p_adams_moulton{i}] = adams_moulton(pd,pdd,p0,t0,dt(i),tmax);
    [t_adams_moulton_linear1{i}, p_adams_moulton_linear1{i}] = adams_moulton_linear1(p0,t0,dt(i),tmax);
    [t_adams_moulton_linear2{i}, p_adams_moulton_linear2{i}] = adams_moulton_linear2(p0,t0,dt(i),tmax);

    % Calculating the exact values.
    tt{i} = [t0:dt(i):tmax];
    pp{i} = 200./(20-10*exp(-7*tt{i}));

    % Finding the real error.
    E_euler(2,i) = getError(p_euler{i},pp{i},dt(i));
    E_heun(2,i)  = getError(p_heun{i},pp{i},dt(i));
    E_implicit_euler(2,i) = getError(p_implicit_euler{i},pp{i},dt(i));
    E_adams_moulton(2,i) = getError(p_adams_moulton{i},pp{i},dt(i));
    E_adams_moulton_linear1(2,i) = getError(p_adams_moulton_linear1{i},pp{i},dt(i));
    E_adams_moulton_linear2(2,i) = getError(p_adams_moulton_linear2{i},pp{i},dt(i));
    
    % Comparison of errors.
    if  i>1
       E_euler(3,i) = E_euler(2,i)/ E_euler(2,i-1);
       E_heun(3,i)  = E_heun(2,i)/ E_heun(2,i-1);
       E_implicit_euler(3,i)  = E_implicit_euler(2,i)/ E_implicit_euler(2,i-1);
       E_adams_moulton(3,i)  = E_adams_moulton(2,i)/ E_adams_moulton(2,i-1);
       E_adams_moulton_linear1(3,i)  = E_adams_moulton_linear1(2,i)/ E_adams_moulton_linear1(2,i-1);
       E_adams_moulton_linear2(3,i)  = E_adams_moulton_linear2(2,i)/ E_adams_moulton_linear2(2,i-1);
    end

    % Finding error based on best approximation.
    % Prune the 'best' results y_*{1} by taking every nth value
    % so we can compare directly against y_*{i}

%     scale = dt(i)/dt(1);
%     E_euler(4,i) = getError(p_euler{i}, p_euler{1}(1:scale:end), dt(i));
%     E_heun(4,i)  = getError(p_heun{i},  p_heun{1}(1:scale:end),  dt(i));
%     E_implicit_euler(4,i) = getError(p_implicit_euler{i}, p_implicit_euler{1}(1:scale:end), dt(i));
%     E_adams_moulton(4,i) = getError(p_adams_moulton{i}, p_adams_moulton{1}(1:scale:end), dt(i));
%     E_adams_moulton_linear1(4,i) = getError(p_adams_moulton_linear1{i}, p_adams_moulton_linear1{1}(1:scale:end), dt(i));
%     E_adams_moulton_linear2(4,i) = getError(p_adams_moulton_linear2{i}, p_adams_moulton_linear2{1}(1:scale:end), dt(i));

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
%disp('Fourth row: The approximation error using the best approximated solution as reference.')

disp(' ')
disp(' ')
fprintf('<strong>Results </strong>')
disp(' ')

E_euler
E_heun
E_implicit_euler
E_adams_moulton
E_adams_moulton_linear1
E_adams_moulton_linear2


method_names = {'Explicit Euler', 'Heun', 'Implicit Euler', 'Adams-Moulton', 'Linearized Adams-Moulton (1)', 'Linearized Adams-Moulton (2)'};
method_names_short = {'ExpEuler', 'Heun', 'ImpEuler', 'AdamsMoulton', 'AdamsMltL1', 'AdamsMltL2'};
method_results = [p_euler; p_heun; p_implicit_euler; p_adams_moulton; p_adams_moulton_linear1; p_adams_moulton_linear2];

% Create plots for each method
plot_method(dt, tt, pp{1}, method_results, method_names, [t0,tmax,0,20], 'kgmbk');

% Print stability table for different criteria
disp 'Bounded error stability criterion';
stability_table(@bounded_error_stability, dt, tt, pp, method_results, method_names_short)
disp 'Falling error stability criterion';
stability_table(@falling_error_stability, dt, tt, pp, method_results, method_names_short)
disp 'Final destination stability criterion';
stability_table(@final_destination_stability, dt, tt, pp, method_results, method_names_short)

