function y = f(x)

y = zeros(size(x));

y(1) = sum(x)    - 3;
y(2) = sum(x.^2) - 7;
y(3) = sum(x.^3) - 11;