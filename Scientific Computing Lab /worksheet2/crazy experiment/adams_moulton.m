
function [t,y]=adams_moulton(f,fd,y0,t0,dt,tmax)
% ADAMS_MOULTON numerically integrates an ODE using the second-
% order Adams-Moulton method.
%   f : Function to be solved
%   y0: Starting value at t0
%   t0: Starting time
%   dt: Timestep size
%   tmax: End time
%
% See also IMPLICIT_EULER, EXPLICIT_EULER, RUNGEKUTTA

t = t0:dt:tmax;
y = zeros(1,length(t));
y(1) = y0;
lin1 = y0;
lin2 = y0;
err = 1e-7;

for i=1:(length(t)-1)
    
    g = @(x) x - y(i) - dt/2 * (f(t(i),y(i)) + f(t(i+1),x));
    gd = @(x) 1 - dt/2 * fd(t(i+1), x);
    
    [y(i+1), lin1, lin2] = newtons_method(g, gd, y(i), err, dt, lin1, lin2);
    lin1
    lin2
        if isnan(y(i+1)) 
        
            sprintf('Newton''s method failed!!!');
            
        end
    
end

