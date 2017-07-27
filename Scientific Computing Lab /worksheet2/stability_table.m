

function result = stability_table(pred, dt, t, exact, approx, methods)
% STABILITY_TABLE generates a table depicting the stability of various
% integration schemes at different timesteps. 
% Letting m := number of methods; n := number of timesteps,
% 
%   pred:    A predicate describing a stability criterion.
%            pred = pred(time_vector, exact_solution_vector, approx_solution_vector)
%
%   dt:      1xn array of timesteps, e.g. [1/32, 1/16, 1/8, 1/2]
%   t:       1xn cell array containing time vectors
%   exact:   1xn cell array containing exact solution vectors
%   approx:  mxn cell array containing approximate solution vectors
%   methods: 1xm cell array of strings naming each method. These need to be
%            valid MATLAB identifiers.
%       
% Returns a table as specified in Worksheet 2.


% Map the matrix of result vectors over the stability predicate
stabilities = zeros([length(methods), length(dt)]);
for method = 1:length(methods)
    for timestep=1:length(dt)
        stabilities(method, timestep) = pred(t{timestep}, exact{timestep}, approx{method, timestep}, false);
    end
end

% Generate headers for each timestep
timestep_headers = cell(size(dt));
for timestep=1:length(dt)
   timestep_headers{timestep} = strcat('dt = 1/', num2str(1/dt(timestep)));
end

% Transpose the matrix to be consistent with the worksheet
% Convert the results into a table for readability
result = array2table(flip(stabilities',1));
result.Properties.VariableNames = methods;
result.Properties.RowNames = flip(timestep_headers);

end

