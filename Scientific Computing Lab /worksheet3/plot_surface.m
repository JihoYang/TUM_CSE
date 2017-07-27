
function plot_surface(X,Y,T,description)
% PLOT_METHOD generates a contour plot and a 3D colored surface.
%
%   X: A mxn matrix of x-values, usually generated via meshgrid()
%   Y: A mxn matrix of y-values, usually generated via meshgrid()
%   T: A mxn matrix of computed temperature values at each grid point
%   description: A short string describing the method used to compute T

    figure;
    hold on;
    [nx,ny] = size(T);

    subplot(1,2,1) 

    surf(X,Y,T);
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    colorbar;
    title([description,' with Nx = ', num2str(nx),', Ny = ',num2str(ny)]);

    subplot(1,2,2)
    contourf(X,Y,T);
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    colorbar;

    box on; grid on;

end

