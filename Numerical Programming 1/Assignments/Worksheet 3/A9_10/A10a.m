function [yy_spline,yy_hermite] = A10a(x,y,w,xx)
% evaluate spline and hermite interpolation for xx.
%   y=f(x) are the function values and
%   w=f'(x) are the values of the derivatives

yy_spline = interp1(x,y,xx,'spline');

yy_hermite = zeros(size(xx));

for j=1:numel(x)-1
    in_interval = xx>=x(j) & xx<=x(j+1);
    yy_hermite(in_interval) = ...
        hermiteinterpol(x(j),x(j+1),y(j:j+1),w(j:j+1),xx(in_interval));
end

end