function [x_iter,omega]=classicIteration(name,A,b,lambda,max_iter)
% (relaxed) Jacobi and SOR iteration methods for solving A*x=b
% function [x_iter,omega]=classicIteration(name,A,b,lambda,max_iter)
% Input:
%       name            'Jacobi' or 'SOR'
%       A,b             (n,n) matrix and (n,1) vector for linear system
%       lambda          relaxation parameter
%       max_iter        number of (maximal) iterations
% Output:
%       x_iter          (n,max_iter) matrix with x-iterates
%       omega           (1,max_iter) backward errors

normA = norm(A,inf); n=numel(b);

switch lower(name)
    case 'sor'
        M=tril(A); N=-triu(A,1); 
    case 'jacobi'
        M=spdiags(diag(A),0,n,n); N=-tril(A,-1)-triu(A,1);
    otherwise
        error(['Unknwon iteration with name: ',name]);
end

x_iter=NaN*ones(n,max_iter); omega=NaN*ones(1,max_iter);
x=zeros(n,1);  Nx=x; r=b; d=M\b;

for k=1:max_iter
    Nx_old=Nx; r_old=r;
    x = lambda*(M\(Nx_old)+d)+(1-lambda)*x;
    nn = normA*norm(x,inf);
    if ~isfinite(nn), return; end   
    x_iter(:,k)=x; Nx = N*x;
    r = Nx-Nx_old + (1-lambda)*r_old;
    omega(k) = norm(r,inf)/nn;    
end


end