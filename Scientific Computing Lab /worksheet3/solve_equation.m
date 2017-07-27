%-----------------------------------------------%
% An implementation of the gauss seidel solver  %
% for our finite difference problem             %
%-----------------------------------------------%

% T     : The solution matrix with boundaries around.
% initial  : Scalar initial guess for each mesh point
% boundary : Sparse matrix with populated edges. size(boundary) = [nx+2,ny+2].
%   - Want to be able to have any boundary value.
% nx     : # of meshpoints in x direction.
% ny     : # of meshpoints in y direction.
% err    : The accepted residual norm.
% f      : Driving function (right hand side of equation).
%               - Matrix or function which takes in f(y,x)
% hx, hy : stepsize in x and y direction.

function T = solve_equation(initial, boundary,err, f,hx,hy)

% n        : # of steps
% nx,ny    : Size of mesh in x and y direction
% R_local  : Local residual
% R_global : Global residual norm

n = 0;
[ny,nx] = size(boundary);
nx = nx-2;
ny = ny-2;
T = full(boundary);
T(2:ny+1,2:nx+1) = initial*ones(ny,nx);
R_global = 2*err;

%dummy variables for cleaner code.
wx = 1/hx^2;  %weight in x direction
wy = 1/hy^2;  %weight in y direction
wh = 2*(1/hx^2 + 1/hy^2); %weight at current location

while R_global>err

    %r: row
    %c: column
    for r = 2:ny+1
        for c = 2:nx+1
            
            T(r,c)= -f(r,c)/wh +...
            wx*(T(r,c-1)+T(r,c+1))/wh +...
            wy*(T(r-1,c)+T(r+1,c))/wh;
        end%for
    end%for

    R_global = get_residual_norm(f,T,wx,wy,wh);
    n = n+1;

end%while

end%function



