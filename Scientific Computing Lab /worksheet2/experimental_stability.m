%-----------------------------------------------------------%
%-- Function that checks if a number solution to a solver --%
%-- is stable based on a very simple criterion            --%
%-----------------------------------------------------------%

function [is_stable, local_error, stable_local_error] = experimental_stability(t, pp, p, visual)
% EXPERIMENTAL_STABILITY is a stability criterion predicate which explores
% the idea of bounding the error using the Jacobian at each timestep. This 
% does not form a useful criterion at the moment, but there may be a good idea 
% buried in here somewhere.

%   t:   Vector of n points in time
%   pp:  Vector of n points of the exact solution
%   p:   Vector of n points of the approximate solution
%   visual:  When true, creates a plot
%

pdd  = @(t,p) 7.*(1-p/5);               

local_error = zeros(size(pp));          % The local error
stable_local_error = zeros(size(pp));   % Tunnel for which the local error has to stay within for us to call it stable
curvature = pdd(t, pp);                 % Actual curvature for every i
scale = 2.1;                            % Scale
stable_minimum_error = 1e-4;
for i=1:length(pp)-1
    local_error(i) = abs(pp(i) - p(i));
    dt = (t(i+1)-t(i));

    % Optimally stable_local_error should be the same kind of function as the local error. 

    % Maybe assume the error is ditributed as some distribution. 
    % From the plots the local error looks something like a boltzmann-maxwell distribution. 
    % This is a bit to much backwards engineering. Will only work for this kind of 
    % function. Ways to make it more general? Make it based on curvature parhaps.

    % This distribution is not the best fit and it should also include the curvature to be applicable
    % to a general function, but is shows the principle.
    stable_local_error(i) = sqrt(2/pi)* t(i)^2 * exp(-t(i)^2/(2*dt^2))/(dt^2) + 1e-4;
end
%-----------------------------% 
%-- Scaling of distribution --%
%-----------------------------% 

% This scaling automaticly accepts the local error
% in the point right after the point which has the largest curvature.

% The error should be biggest around the largest curvature.
[maximum, index] = max(abs(curvature));

% The largest error will be shown in the next point.
index = index + 1; 

% Scale the threshold. Also multiplying by 1.2 to make it slightly bigger than the local error
% at index.
stable_local_error = 1.2*stable_local_error * local_error(index) /stable_local_error(index)
%----------------------------------------------------%

stability = local_error < stable_local_error;
if (visual)
    figure;
    hold on;
    semilogx(t, local_error,'b');
    semilogx(t, stable_local_error,'r');
    semilogx(t, stability,'y');
    legend('local error', 'threshold', 'stability (0: unstable, 1:stable)')
end
is_stable = all(stability);
end
