


function [is_stable, transient_mask, err_delta] = falling_error_stability(t, y_exact, y_appx, visual)
% FALLING_ERROR_STABILITY is a stability criterion predicate. The basic idea is that as 
% long as our y_appx is within some neighborhood of a stable critical point,
% the stepwise global error should constantly decrease. This does not always 
% hold true for higher-order methods with smaller timesteps in transient regions,
% so we ignore the transient region just like we did in BOUNDED_ERR_STABILITY.
%
% t:       Array of time points [t_0, t_1, ..., t_max]
% y_exact: Array of exact y
% y_appx:  Array of approximate y
% visual:  When true, creates a plot
% 
% Returns 0 or 1 depending on stability, along with debugging information
% See also: BOUNDED_ERROR_STABILITY, FINAL_DESTINATION_STABILITY

pd = @(t,p) 7.*(1-p/10).*p;

% Create a mask to identify 'transient' regions := points where dy/dt>threshold
transient_mask = abs(pd(t,y_exact)) > 0.5;
% We use y_exact here because otherwise, when y_approx is very bad, we would 
% mistakenly identify the entire curve as being 'transient'

% Error delta expresses |y_exact(n) - y_approx(n)| - |y_exact(n+1) - y_approx(n+1)|
err = abs(y_exact - y_appx);
err_delta = err(2:end-1) - err(3:end);
falling_err_mask = err_delta > -1e-6;

% We always need to skip the first point since y_exact(1) = y_approx(1)
transient_mask = transient_mask(3:end);

% Stability => for all points p, (error constantly decreases) OR (p is transient)
stable_mask = falling_err_mask | transient_mask;
is_stable = all(stable_mask); 

if(visual)
    figure;
    subplot(2,1,1);
    hold on;
    title('Actual y');
    plot(t, y_exact, '-');
    plot(t(~stable_mask), y_appx(~stable_mask), 'r*');
    plot(t(falling_err_mask), y_appx(falling_err_mask), 'ko');
    plot(t(transient_mask), y_appx(transient_mask), 'k+');
    subplot(2,1,2);
    hold on;
    title('Absolute error in y');
    plot(t, err, '-');
    plot(t(~stable_mask), err(~stable_mask), 'r*');
    plot(t(falling_err_mask), err(falling_err_mask), 'ko');
    plot(t(transient_mask), err(transient_mask), 'k+');
end

end

