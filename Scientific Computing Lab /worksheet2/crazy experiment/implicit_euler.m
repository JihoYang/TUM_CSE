%---------------------------------------------%
%-- Implicit Euler ODE Approximation Solver --%
%---------------------------------------------%
%---------Solves the ODE: fd = f(t,y)---------%
%---------------------------------------------%

% f :   Function to be solved
% fd:   Derivatve of f
% y0:   Starting value at t0
% t0:   Starting time
% dt:   Timestep size
% tmax: Endtime
% t:    Timesteps
% err:  Accuracy limit

function [t,y]=implicit_euler(f,fd,y0,t0,dt,tmax)

%Initialise time and corresponding y (output) matrices to which the
%computed solutions will be stored

t = [t0:dt:tmax];
y = zeros(1,length(t));

%Set initial condition and accuracy limit (for Newton's method)

y(1) = y0;
err = 1e-7;

%Loop for implicit Euler method. The method implicitely solves the solution
%at the next time step (y(i+1)) by using the initial value (or solution at
%current time step y(i)) and the gradient at the next time step
%(f(t+1,y+1)). In order to compute y(i+1) (which is used to compute
%f(t+1,y+1)), Newton's method is implemented. A new function g and
%derivative of g are defined here and used as input function handles for
%Newton's method 

    for i=1:(length(t)-1)
        
        g=@(y_new) y_new-y(i)-dt*f(t(i+1),y_new);
        gd=@(y_new) 1-dt*fd(t(i+1),y_new);
        [y(i+1), lin1, lin2] =newtons_method(g,gd,y(i),err,dt, y0, y0);
        
        if isnan(y(i+1)) 
        sprintf('Newton''s method failed!!!');
        
    end
    
end

