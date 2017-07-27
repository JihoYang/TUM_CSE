
function T = gauss_seidel(T_old, initial, boundary,err, f,dt,hx,hy)
% GAUSS_SEIDEL is a Gauss-Seidel solver customized for the equation
% (I-dt*L) T = T_old + dt*f, which comes from the implicit 
% Euler formulation of the heat equation. 
%
% 
% T_old : 
% initial  : Scalar initial guess for each mesh point
% boundary : Sparse matrix with populated edges. size(boundary) = [nx+2,ny+2].
%   - Want to be able to have any boundary value.
% nx     : # of meshpoints in x direction.
% ny     : # of meshpoints in y direction.
% err    : The accepted residual norm.
% f      : Driving function (right hand side of equation).
%               - Matrix or function which takes in f(y,x)
% hx, hy : stepsize in x and y direction.
% T     : The solution matrix with boundaries around.

n = 0;
[ny,nx] = size(boundary);
nx = nx-2;
ny = ny-2;
T = full(boundary);
T(2:ny+1,2:nx+1) = initial*ones(ny,nx);
R_global = 2*err;

wx = -dt/(hx^2);  %weight in x direction
wy = -dt/(hy^2);  %weight in y direction
wh = 1+(2*dt*(1/hx^2 + 1/hy^2)); %weight at current location
b = dt*f + T_old;

while R_global>err

    for r = 2:ny+1 %r: row
        for c = 2:nx+1 %c: column
            
            T(r,c)=  (wx * (T(r,c-1) + T(r,c+1)) ...
                    + wy * (T(r-1,c) + T(r+1,c)) ...
                    - b(r,c)) / -wh;
        end
    end

    R_global = get_residual_norm(b,T,wx,wy,wh);
    n = n+1;

end

end



