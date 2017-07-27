

function result = implicit_euler_gauss_seidel(t0, dt, tmax, initial, boundary, forcing)
% IMPLICIT_EULER_GAUSS_SEIDEL uses an implicit Euler scheme with a custom Gauss-Seidel
% solver to integrate the PDE du/dt = L(u)+f, where L is the Laplacian. 
% This assumes a uniform rectangular mesh and Dirichlet boundary conditions.
% 
%   t0:        The starting time
%   dt:        The timestep size
%   tmax:      The stopping time
%   initial:   a (m+2)*(n+2) matrix containing initial conditions at each grid point
%   boundary:  a (m+2)*(n+2) sparse matrix containing Dirichlet boundary 
%              conditions at the edges of the grid.
%   forcing:   a (m+2*n+2) matrix containing the forcing term at each grid point
%
%
%   returns:   a struct containing a (m+2)*(n+2) matrix representing the 
%              solution (including boundary conditions) 


[ny,nx] = size(forcing);
nx = nx-2;
ny = ny-2;
hx = 1/(nx+1);
hy = 1/(ny+1);
solution = initial;
for t = t0+dt:dt:tmax
    
    % solution = (I-dt*L) \ (solution + dt*forcing);
    % Overwrite solution on each iteration to save memory
    solution = gauss_seidel(solution,1,boundary,1e-6,forcing,dt,hx,hy);
end

result.method = 'Implicit Euler Gauss-Seidel';
result.n = size(initial,1);
result.dt = dt;
result.t = tmax;
result.solution = solution;
result.eigMin = 0;
result.eigMax = 0;