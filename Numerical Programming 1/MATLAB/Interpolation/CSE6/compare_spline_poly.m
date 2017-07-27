function compare_spline_poly

% (c) Folkmar Bornemann 26.1.2005

% parameter

n_max = 200;        % maximal polynomial degree 

% Runge function (Runge 1901: alpha = 5)

alpha = 5;
f  = @(x) 1./(1+alpha^2*x.^2);

% theoretical prediction of convergence

s_cheb = sqrt(1+1/alpha^2); % abscissa
K = s_cheb+1/alpha;         % convergence constant  

% plots of interpolants

fig = initdraw;
for n=1:n_max
    [p,E1(n)] = DrawInterp(f,-1,1,n,'spline',1);
    [p,E2(n)] = DrawInterp(f,-1,1,n,'lagrange',2);
end

% error plots

ErrorPlot(E1,4,'spline',99,1)
ErrorPlot(E2,K,'lagrange',99,2)

% =========================================================================

function [p,err] = DrawInterp(f,a,b,n,method,sub)

% parameter --------------------------------------------------------

N = 2000;                       % number of points for plotting/error
n_max = 100;                    % maximum degree to be plotted
y_min = -0.5; y_max = 1.5;      % y range
x_str = -0.2; y_str = 1.125;    % position of text

% ------------------------------------------------------------------

x = linspace(a,b,N);

switch method
    case 'spline'   
        knots = linspace(a,b,n+1);
        fval = f(knots);
        pp =spline(knots,fval);
        title_str='not-a-node cubic spline interpolation (equidistant)';
        p = ppval(pp,x);
    case 'lagrange' 
        knots = (a+b)/2+(b-a)/2*cos(linspace(0,pi,n+1));
        fval = f(knots);
        lambda = weights(knots); 
        p = barycentric2(x,knots,fval,lambda);
        title_str='barycentric formula (Chebyshev nodes)';
end

if n <= n_max
    subplot(1,2,sub);
    h = plot(x,p,'b-',knots,fval,'r.'); 
    set(h,'LineWidth',2,'MarkerSize',16);
    str = sprintf('n = %i',n);
    h=text(x_str,y_str,str);
    set(h,'Fontsize',[24],'Color','black');
    xlabel('x');
    ylabel('f(x)');
    title(title_str);
    axis([a,b,y_min,y_max]); drawnow;
end

err = norm(f(x)-p,inf);

% =========================================================================

function ErrorPlot(E,K,method,fig,sub)

h=figure(fig); set(h,'Name','error plot')
subplot(1,2,sub);

m = length(E);
switch method
    case 'lagrange'
        h=semilogy(2:m+1,E,'b-'); hold on
        set(h,'LineWidth',2,'MarkerSize',16);
        p = polyfit(1:20,log(E(1:20)),1);
        C= exp(p(2));
        h=semilogy(1:m,C./K.^(1:m),'r-');
        set(h,'LineWidth',2,'MarkerSize',16);
        axis([1 m+1, 10^-15, 10^0])
        title('barycentric formula (Chebyshev nodes)')
        xlabel('dimension');
        ylabel('absolute error');
    case 'spline'
        h = loglog(2:m+1,E,'b-'); hold on
        set(h,'LineWidth',2,'MarkerSize',16);
        p = polyfit(log(125:150),log(E(125:150)),1);
        C= exp(p(2));
        h=semilogy(1:m,C./(1:m).^K,'r-');
        set(h,'LineWidth',2,'MarkerSize',16);
        axis([1 m+1, 10^-15, 10^0])
        title('cubic spline interpolation (equidistant)')
        xlabel('dimension');
        ylabel('absolute error');
end

% =========================================================================

function fig = initdraw

close all
screensize=get(0,'ScreenSize');
fig=figure('Outerposition',screensize,'Units',get(0,'Units'));
figure(fig);
set(fig,'DoubleBuffer','on','Name','Splines vs. Polynome');

