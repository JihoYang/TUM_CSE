
%--------------------------%
%   Scientific Computing   %
%        Worksheet 4       %
%--------------------------%
%        Nathan Brei       %
%         Jiho Yang        %
%   Tord Kriznik Sørensen %
%--------------------------%

close all;

% Create a test matrix indexed by time and space discretizations

t0 = 0;
tmax = 10;

methods = {@explicit_euler_direct, @adams_moulton, @implicit_euler_direct};
methods_s = [ 'Adams-Moulton '; 'Explicit Euler'; 'ImpliciEulerGS'];
timesizes = (1/2).^(6:12);
spacesizes = [3,7,15,31];
plot_times = [0,1] / 8;

test_matrix = cell((length(plot_times)-1)*length(timesizes)*length(spacesizes)*length(methods), 10);

% Allocate cell arrays for problem inputs
[mesh, boundary, forcing, initial] = deal(cell(length(spacesizes), 1));

% Set up the problem inputs for each space discretization
for s=1:length(spacesizes)

    % Create rectangular grid
    n       = spacesizes(s);
    h       = 1/(n+1);
    [X,Y]   = meshgrid(0:h:(n+1)*h);
    mesh{s} = {X;Y};

    % Specify forcing function (RHS) along grid
    forcing{s} = zeros(n+2);

    % Specify Dirichlet boundary conditions along outside of grid
    boundary{s}        = sparse(n+2);
    boundary{s}(:,1)   = 0;
    boundary{s}(:,n+2) = 0;
    boundary{s}(1,:)   = 0;
    boundary{s}(n+2,:) = 0;

    % Specify initial solution T(t=0)
    initial{s} = zeros(n+2);
    initial{s}(2:n+1, 2:n+1) = ones(n);
end

% Initiating variables for plotting.
counter = 0; % Location to plot next subfigure
lm      = length(methods);
lpt     = length(plot_times)-1;

% Populate test matrix with results for each (spacesize, timesize) pair
output_row=1;
tic;
for m = 1:length(methods)
    for s = 1:length(spacesizes)
        for t = 1:length(timesizes)

            last_solution = initial{s};

            for p = 2:length(plot_times)

                result = methods{m}(plot_times(p-1), ...
                                    timesizes(t),    ...
                                    plot_times(p),   ...
                                    last_solution,   ...
                                    boundary{s},     ...
                                    forcing{s});
                
                % Calculate error metrics
                result.energyChange = energy_change(last_solution, result.solution);
                result.max = max(max(result.solution));
                result.stable = all(all(result.solution<=1 & result.solution>=0));

                % Use solution from last interval as initial condition for next
                last_solution = result.solution;
                
                test_matrix(output_row,:) = struct2cell(result);
                output_row = output_row+1;
            end            
        end
    end

end
time_solve = toc;

results = cell2table(test_matrix, 'VariableNames', fieldnames(result));

% --------------- %
%     Plotting    % 
% --------------- %
% Sort results for for plotting
results = sortrows(results,'t');
results = sortrows(results, 'method')


old = 0;
ls = length(spacesizes);
lt = length(timesizes);

n = results.n((old+1):(old+ls*lt));    % Same for all timeslices
t_length = (old+ls*lt);                % Same for all timeslices
dt  = results.dt((old+1):(old+ls*lt)); % Same for all timeslices

tic;
for m = 1:lm
    for p = 1:lpt
        fig = p + lpt*(m-1);
        sub_title = ['t = ', num2str(plot_times(p+1)) , ', method: ', ... 
                    methods_s(m,:)];
        T  = results.solution((old+1):(old+ls*lt)); 

        plot_timeslice(mesh,        ...
                         n,           ... 
                         t_length,    ...
                         T,           ...
                         ls, lt,      ...
                         fig,         ...
                         sub_title);

        old = old+ls*lt;

    end


end
time_plot = toc;
