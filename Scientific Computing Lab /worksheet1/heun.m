function [t,y] = heun(f, y0, t0, dt, tmax)
% HEUN computes 
%   f : Function to be solved
%   y0: Starting value at t0
%   t0: Starting time
%   dt: Timestep size
%   tmax: End time
%
% See also EXPLICIT_EULER, RUNGEKUTTA

% Allocate vectors 
t = t0:dt:tmax;
y = zeros(size(t));
y(1) = y0;

% We separate iterating over t from the scheme itself for conceptual clarity.
% The compiler will inline the function so there is no performance cost.
for i = 1:length(t)-1
    y(i+1) = heun_step(f, y(i), t(i), dt);     
end

end


function y = heun_step(f, x, t, dt)
% HEUN_STEP is the Heun method as described in the 
% Scientific Computing lecture slides.

ydot = f(t + dt, x + dt * f(t, x));
y = x + dt * (1/2) * (f(t,x) + ydot);

end
