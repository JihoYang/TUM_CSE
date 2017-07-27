function result = gaussSeidel(forcing, boundary)
% GAUSSSEIDEL solves the discrete Poisson equation on a rectangular grid. 
% It uses a Gauss-Seidel iterative solver in order to minimize memory use 
% on large grids.
% 
%   forcing:   a m*n matrix containing the forcing term at each grid point
%   boundary:  a (m+2)*(n+2) sparse matrix containing Dirichlet boundary 
%              conditions at the edges of the grid.
%
%   returns:   a struct containing a (m+2)*(n+2) matrix representing the 
%              solution (including boundary conditions) along with
%              metadata, including the runtime and the memory required to
%              represent the Laplacian operator.


    tic;

    err      = 1e-4;
    initial  = 0;

    [ny,nx] = size(forcing);
    hx = 1/(nx+1);
    hy = 1/(ny+1);
    
    % solve_equation expects a forcing matrix that includes boundaries
    forcing_full = zeros(size(forcing)+2);    
    forcing_full(2:ny+1,2:nx+1) = forcing;
    
    solution = solve_equation(initial,boundary,err,forcing_full,hx,hy);
    solution = solution(2:ny+1, 2:nx+1);
    runtime  = toc;

    result          = struct();
    result.method   = 'Gauss-Seidel';
    result.size     = nx;
    result.solution = solution;
    result.runtime  = runtime;
    result.storage  = numel(forcing_full);

end
