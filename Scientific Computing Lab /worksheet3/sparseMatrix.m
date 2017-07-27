function result = sparseMatrix(forcing, boundary)
% SPARSEMATRIX solves the discrete Poisson equation on a rectangular grid. 
% It uses a sparse matrix representation of the Laplacian and MATLAB's
% backslash solver.
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

    [nx, ny] = size(forcing);
    
    A = make_laplacian_matrix_sparse(nx,ny);
    
    % Reshape forcing function from matrix to column vector
    b = reshape(forcing', [numel(forcing') 1]);
    
    % TODO: sparseMatrix assumes boundary = 0 for now

    t = A\b;
    
    % Reshape solution temperatures from column vector to matrix
    T = reshape(t, [nx ny]);
    
    runtime = toc;
    
    result          = struct();
    result.method   = 'Sparse matrix';
    result.size     = nx;
    result.solution = T;
    result.runtime  = runtime;
    result.storage  = nnz(A) + numel(t) + numel(b);

end
