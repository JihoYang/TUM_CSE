function [mu,v,omega,k] = PowerMethod(A,v,tol)

[m,~] = size(A);            % dimension
normA = norm(A,'fro');      % matrix norm
omega = inf;                % initialize backward error				
w = A*v;                    % first A*v
k = 0;
while omega > tol           % backward error too large?
    k = k + 1;              % step count
	v = w/norm(w);          % normalize A*v			
	w = A*v;                % new A*v	
	mu = v'*w;              % Rayleigh quotient		
	r = w - mu*v;           % residual               
	omega = norm(r)/normA;  % backward error
end
