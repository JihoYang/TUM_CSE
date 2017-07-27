rho = 10/n;

rng(843);
A = sprandn(n,n,rho) + 4*speye(n);
b = ones(n,1);