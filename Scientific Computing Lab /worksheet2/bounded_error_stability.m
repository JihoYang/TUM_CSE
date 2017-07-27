

function [is_stable, transient_mask, bounded_mask] = bounded_error_stability(t, y_exact, y_appx, visual)
% BOUNDED_ERROR_STABILITY is a stability criterion predicate. The basic idea 
% is that as long as our y_appx is within some neighborhood of a stable critical point,
% the (stepwise) global error should remain bounded. We establish two
% thresholds in order to define a concept of 'close to zero' for f(t,p) and 
% |y_exact-y_appx|. Stability implies that when f->0, the error->0.
%
% t:       Array of time points [t_0, t_1, ..., t_max]
% y_exact: Array of exact y
% y_appx:  Array of approximate y
% visual:  When true, creates a plot
% 
% Returns 0 or 1 depending on stability, along with debugging information
% See also: FALLING_ERROR_STABILITY, FINAL_DESTINATION_STABILITY


pd = @(t,p) 7.*(1-p/10).*p;

% Create a mask to identify 'transient' regions := points where dy/dt>threshold
transient_mask = abs(pd(t,y_exact)) > 0.5;
% We use y_exact here because otherwise, when y_approx is very bad, we would 
% mistakenly identify the entire curve as being 'transient'

% Create another mask to identify regions where error is within threshold
err = abs(y_exact-y_appx);
bounded_mask = err < 1;

% Stability => for all points p, (error is bounded at p) OR (p is transient)
stable_mask = bounded_mask | transient_mask;
is_stable = all(stable_mask);

if(visual)
    figure;
    subplot(2,1,1);
    hold on;
    title('Actual y');
    plot(t, y_exact, '-');
    plot(t(~stable_mask), y_appx(~stable_mask), 'r*');
    plot(t(bounded_mask), y_appx(bounded_mask), 'ko');
    plot(t(transient_mask), y_appx(transient_mask), 'k+');
    subplot(2,1,2);
    hold on;
    title('Absolute error in y');
    plot(t, err, '-');
    plot(t(~stable_mask), err(~stable_mask), 'r*');
    plot(t(bounded_mask), err(bounded_mask), 'ko');
    plot(t(transient_mask), err(transient_mask), 'k+');
end

end

