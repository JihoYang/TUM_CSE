%------------------------------------------%
% Implementation of the error function
% E = sqrt(dt/5* sum(p_k - p_(k,exact)))
%------------------------------------------%


% f: approximate solution from a solver
% ff: exact or best solution available
% dt: timestep used to get f

function E = getError(f,ff,dt)
    E = sqrt(dt/5 * sum((f-ff).^2));    
end

