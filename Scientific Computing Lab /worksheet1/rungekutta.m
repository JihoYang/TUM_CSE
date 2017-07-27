
%------------------------------------%
%------------------------------------%
% Runge-Kutta 4th order (RK4) solver %
%------------------------------------%
%------------------------------------%
%---- Solves the ODE: fd = f(t,y)----%
%------------------------------------%
%------------------------------------%


% y0:   Starting value at t0
% t0:   Starting time
% dt:   Timestep size
% tmax: Endtime
% t:    Timesteps

% Y1: Slope at the starting point.
% Y2: Slope at the middle point. The middle point is found by going along the slope Y1 from yn until
%     half a timestep is reached. 
% Y3: Slope at the middle point. The middle point found as for Y2, but by using Y2 instead of Y1.
% Y4: Slope at the end point. The end point is found by going along the slope Y3 from yn until
%     the next timestep is reached. 
% y:  The vector of points found by the method. Each point is found by going from the previous
%     point to the next along a slope that is a weighted sum of Y1, Y2, Y3 and Y4.

% The expression for y(i+1) is found by expanding it into a taylor series with
% five terms and then assume that y(i+1) = y(i) + (a1*Y1 + a2*Y2 + a3*Y3 +a4*Y4)*dt.
% By equating the two terms for y(i+1) we get the scheme for RK4.


function [t, y] = rungekutta(f,y0,t0,dt,tmax)

    t = [t0:dt:tmax];
    y = zeros(1,length(t));
    y(1) = y0;
    
    for i=1:(length(t)-1)

        Y1 = f(t(i),y(i));
        Y2 = f(t(i)+ 0.5*dt,y(i) + 0.5*dt*Y1);
        Y3 = f(t(i)+ 0.5*dt,y(i) + 0.5*dt*Y2);
        Y4 = f(t(i) + dt,y(i) + dt*Y3);
        
        y(i+1) = y(i) + dt/6 * (Y1 + 2*Y2 + 2*Y3 +Y4);

     end
end
