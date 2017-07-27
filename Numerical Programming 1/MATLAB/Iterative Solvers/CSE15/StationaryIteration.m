tic, 

% Gauss-Seidel:

% N = -triu(A,1);
% M = spdiags(diag(A),0,n,n) + tril(A,-1);

% Jacobi

N = -tril(A,-1)-triu(A,1);
M = spdiags(diag(A),0,n,n);


K = 10000; tol = 1e-10;

d = M\b; 
x = zeros(n,1); normA = norm(A,inf);
for k = 1:K
    x = M\(N*x) + d;
    r = b - A*x;
    omega = norm(r,inf)/(normA*norm(x,inf));
    if omega < tol, break; end
end
toc

[k,omega]
