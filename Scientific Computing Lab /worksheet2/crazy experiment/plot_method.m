
function plot_method(dt, t, exact, appx, methods, bounds, colors)
% PLOT_METHOD generates a plot comparing the exact solution to an IVP
% against approximations calculated at different timesteps. 
% Letting m := number of methods; n := number of timesteps,
% 
%   dt:      1xn array of timesteps, e.g. [1/32, 1/16, 1/8, 1/2]
%   t:       1xn cell array of time vectors
%   exact:   Exact solution vector, corresponding to t{1}
%   approx:  mxn cell array containing approximate solution vectors
%   methods: 1xm cell array of strings naming each method. 
%   bounds:  The axis bounds [x_min,x_max,y_min,y_max]
%   colors:  1xn cell array (or string) of color designators

for method=1:length(methods)

    figure;
    hold on;
    title(strcat(methods(method), ' Method for Different Timesteps (dt)'));
    legend_entries = cell(1, length(dt)+1);

    plot(t{1}, exact, 'r', 'LineWidth', 1.5);
    legend_entries{1} = 'Analytical Solution';

    for timestep=1:length(dt)
        plot(t{timestep}, appx{method,timestep},'o--','LineWidth',1.5,'Color',colors(timestep));
        legend_entries(timestep+1) = strcat(methods(method),' (dt=1/', num2str(1/dt(timestep)),')');
    end

    set(gca, 'FontSize', 15);
    axis(bounds);
    xlabel('t','FontSize',20);
    ylabel('P(t)','FontSize',20);
    box on; grid on;
    legend(legend_entries, 'Location', 'northeast');

end

