
%---------------------------------------------%
%-- Explicit Euler ODE Approximation solver --%
%---------------------------------------------%
%---------Solves the ODE: fd = f(t,y)---------%
%---------------------------------------------%

% f :   Function to be solved
% y0:   Starting value at t0
% t0:   Starting time
% dt:   Timestep size
% tmax: Endtime
% t:    Timesteps

function [t,y]=explicit_euler(f,y0,t0,dt,tmax)

%Initialise time and corresponding y (output) matrices to which the
%computed solutions will be stored

t = [t0:dt:tmax];
y = zeros(1,length(t));
y(1) = y0;

%Loop for Explicit Euler approximation. The method explicitely uses the
%solution at current time step (y(i)) and the corresponding gradient
%f(t(i),y(i)) to compute the solution (y(i+1)) at the next time step
%(t0+dt)

    for i=1:(length(t)-1)
        
        y(i+1)=y(i)+dt*f(t(i),y(i));
        
    end
end

