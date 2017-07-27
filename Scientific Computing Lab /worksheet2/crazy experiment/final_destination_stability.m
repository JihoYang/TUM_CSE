

function is_stable = final_destination_stability(t, y_exact, y_appx, visual)
% FINAL_DESTINATION_STABILITY is a stability criterion predicate. 
% The basic idea is that if our (method, ode, timestep) combination is stable,
% the final y_n will settle down to a value close to that of the analytical 
% solution. This assumes that our ODE has large stable critical points and 
% t_max is large enough for our curve to reach them. It also assumes that 
% global truncation error remains small. 
%
% t:       Array of time points [t_0, t_1, ..., t_max]
% y_exact: Array of exact y
% y_appx:  Array of approximate y
% visual:  When true, creates a plot
% 
% Returns 0 or 1 depending on stability
% See also: BOUNDED_ERROR_STABILITY, FALLING_ERROR_STABILITY


is_stable = abs(y_exact(end)-y_appx(end)) < 0.01;

