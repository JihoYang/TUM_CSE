
function [t,y_linear1]=adams_moulton_linear1(y0,t0,dt,tmax)

% ADAMS_MOULTON_linear1 numerically integrates an ODE using linearised
% Adams-Moulton method.

%   y0: Starting value at t0
%   t0: Starting time
%   dt: Timestep size
%   tmax: End time

% See also ADMAS_MOULTON, IMPLICIT_EULER, EXPLICIT_EULER

t = t0:dt:tmax;
y_linear1 = zeros(1,length(t)); 
y_linear1(1) = y0;

for i=1:(length(t)-1)
    
    y_linear1(i+1) = (1+7/2*dt*(2-y_linear1(i)/10))*y_linear1(i)/(1+7/20*dt*y_linear1(i));
    
end
