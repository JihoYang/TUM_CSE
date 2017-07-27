function [x,fail,steps] = NewtonRelax(f,Df,x,tol)

%		x           solution
%		fail        flag: 
%                       0 = convergence
%                       1 = relaxation too small
%                       2 = too many iterations
%                       3 = singular Jacobian Df
%		steps       number of iterations
			
lambda_min = 1e-6;			                % parameters
max_steps = 125;

fail = 2;				                    % initialization
steps = 0; 
lambda = 1;
y = f(x); J = Df(x);                        % evaluation

while fail==2 && steps<max_steps 	        % main loop
    steps = steps+1; 
    x_old = x; 
    [L,U] = lu(J);                          % LU decomposition of Jacobian
    dx = -U\(L\y); dx_norm = norm(dx);      % Newton correction
    if isnan(dx_norm)                       % singular?
        fail = 3; break
    end
    while true                              % relaxation loop
        x = x_old + lambda*dx; 
        if dx_norm <= tol  		            % convergence?
            fail = 0; break;
        end
        y = f(x); J = Df(x);                % evaluation
        dx_bar = -U\(L\y);                  % simplified Newton correction
        if norm(dx_bar) < dx_norm 	        % natural monotonicity?
            lambda = min(1,2*lambda);       % increase relaxation
            break;
        end
        lambda = lambda/2;                  % decrease relaxation
        if lambda < lambda_min              % too small?
            fail = 1; break;
        end
    end
end



