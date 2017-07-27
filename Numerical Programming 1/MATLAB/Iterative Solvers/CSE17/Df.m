function J = Df(x)

n = length(x);
J = zeros(n,n);

J(1,:) = ones(1,3);
J(2,:) = 2*x.';
J(3,:) = 3*(x.^2).';

