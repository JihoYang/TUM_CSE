
function [t,y_linear2]=adams_moulton_linear2(y0,t0,dt,tmax)

% ADAMS_MOULTON_linear2 numerically integrates an ODE using linearised
% Adams-Moulton method.
%   y0: Starting value at t0
%   t0: Starting time
%   dt: Timestep size
%   tmax: End time

% See also ADMAS_MOULTON, IMPLICIT_EULER, EXPLICIT_EULER

t = t0:dt:tmax;
y_linear2 = zeros(1,length(t));
y_linear2(1) = y0;

for i=1:(length(t)-1)
    
    y_linear2(i+1) = (y_linear2(i)+dt/2*(7*y_linear2(i)*(1-y_linear2(i)/10)))/(1-7/2*dt*(1-y_linear2(i)/10));
    
end