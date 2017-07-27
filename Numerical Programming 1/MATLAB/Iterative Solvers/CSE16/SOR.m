function [x,omega,k] = SOR(A,b,tol,lambda)

[n,n] = size(A);

N = -triu(A,1);
M = spdiags(diag(A),0,n,n) + tril(A,-1);


K = 10000;

d = M\b; 
x = zeros(n,1); normA = norm(A,inf);
for k = 1:K
    x = lambda*(M\(N*x) + d) + (1-lambda)*x;
    r = b - A*x;
    omega = norm(r,inf)/(normA*norm(x,inf));
    if omega < tol, break; end
end
