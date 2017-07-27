
%---------------------------%
% This function returns the %
% sparse laplacian matrix   %
%---------------------------%

% nx: # of entries in x-direction
% ny: # of entries in y-direction

function A = make_laplacian_matrix_sparse(nx,ny)

    hx = 1/(nx+1);
    hy = 1/(ny+1);

    wx = 1/hx^2;
    wy = 1/hy^2;
    wh = -2*(wx + wy);
    
    n = nx*ny;

    D  = sparse(1:n,1:n,wh*ones(1,n),n,n);
    D1 = sparse(2:n,1:n-1,wx*ones(1,n-1),n,n);
    D2 = sparse((nx+1):n,1:n-nx,wy*ones(1,n-nx),n,n);

    for i = 1:(nx-1)
        D1(i*nx+1, i*nx) = 0;
    end

    A =D2 + D1 + D + D1' + D2';

end
