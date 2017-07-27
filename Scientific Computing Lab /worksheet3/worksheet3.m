
%--------------------------%
%   Scientific Computing   %
%        Worksheet 3       %
%--------------------------%
%        Nathan Brei       %
%         Jiho Yang        %
%   Tord Kriznik Sørensen %
%--------------------------%


close all;


% Create a test matrix for a range of methods and granularities
methods = {@sparseMatrix, @fullMatrix, @gaussSeidel};
granularities = [7 15 31 63 ];

test_matrix = cell(length(methods)*length(granularities), 5);

% Allocate cell arrays for problem inputs
[mesh, boundary, forcing, exact_soln] = deal(cell(length(granularities), 1));

% Set up the problem inputs for each step size
for s=1:length(granularities)
    
    % Create rectangular grid
    n = granularities(s);
    h = 1/(n+1);
    [X,Y] = meshgrid(h:h:n*h);
    mesh{s} = {X;Y};
    
    % Discretize forcing function (RHS) along grid
    forcing{s} = -2*(pi^2)*sin(pi*X).*sin(pi*Y);

    % Specify Dirichlet boundary conditions along outside of grid
    boundary{s} = sparse(n+2); 
    boundary{s}(:,1) = 0;
    boundary{s}(:,n+2) = 0;
    boundary{s}(1,:) = 0;
    boundary{s}(n+2,:) = 0;

    % Discretize analytical solution along grid to compare
    exact_soln{s} = sin(pi*X).*sin(pi*Y);
end

% Populate test matrix with results for each (method, size) pair
for m=1:length(methods)
    for s=1:length(granularities)
        result = methods{m}(forcing{s}, boundary{s});
        test_matrix((m-1)*length(granularities)+s,:) = struct2cell(result);
    end
end


%------------------------------------------%
% Display test matrix as a formatted table %
%------------------------------------------%

disp('WARNING: On some systems, this script must be run several times to achieve accurate runtimes.')
disp('')
disp('Runtime is measured in seconds')
disp('Storage is measured in number of doubles')
disp('Storage does not take variables that are used in all solvers into account')

performance = cell2table(test_matrix, 'VariableNames', fieldnames(result))

%----------------------------%
% Plot the analytic solution %
%----------------------------%

plot_surface(mesh{end}{1}, mesh{end}{2}, exact_soln{end}, 'Analytic solution');

%-------------------------------------------------%
% Populate error table for Gauss-Seidel solutions %
%-------------------------------------------------%

gs = performance(strcmp(performance.method, 'Gauss-Seidel'), :);
rows = height(gs);
[error, reduction] = deal(zeros(rows,1));

for row=1:rows
    error(row) = get_error(gs.solution{row},exact_soln{row});

    if row>1
        reduction(row) = error(row-1)/error(row);
    end
    plot_surface(mesh{row}{1}, mesh{row}{2}, gs.solution{row}, gs.method{row});    
end
errors = table(gs.size, error, reduction, 'VariableNames', {'size', 'error', 'reduction'})

tradeoffs


