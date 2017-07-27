function [mu,v,omega,k] = InversePowerMethod(A,muShift,v,tol)

[m,~] = size(A);                            % dimension
normA = norm(A,'fro');                      % matrix norm
omega = inf;                                % initialize backward error
[L,U,p] = lu(A - muShift*eye(m),'vector');  % keep LU decomposition 
k = 0;
while omega > tol					% backward error too large?
    k = k + 1;                      % step count
	w = U\(L\v(p));					% solve linear system
	normW = norm(w);				% norm of new vector
	z = v/normW;					% auxiliary vector
	v = w/normW;					% normalize
	rho = v'*z;						% auxiliary quantity 
	mu = muShift + rho;             % Rayleigh quotient
	r = z - rho*v; 					% residual
	omega = norm(r)/normA;          % backward error
end