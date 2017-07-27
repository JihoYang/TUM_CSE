function RungePhenomenon

% (c) Folkmar Bornemann 26.1.2005

% parameter

n_max = 200;        % maximum degree

% Runge function (Runge 1901: alpha = 5)

alpha = 5;
f  = @(x) 1./(1+alpha^2*x.^2);

% abscissa of convergence domain for equidistant nodes

r1 = @(x) -2+2*x*atan(1/x)+log(1+x^2);
r2 = @(x) -2-(-1+x)*log(1-x)+(1+x)*log(1+x)-r1(1/alpha);
s_equi = abs(fzero(r2,0.5));    % abscissa

% abscissa of convergence domain for Chebyshev nodes

s_cheb = sqrt(1+1/alpha^2); % abscissa
K = s_cheb+1/alpha;         % convergence constant  

% plot of interpolants

fig = initdraw;
for n=1:n_max
               DrawInterp(f,-1,1,n,'equi',s_equi,1);
    [p,E(n)] = DrawInterp(f,-1,1,n,'cheb',s_cheb,2);
end

% error plot for Chebyshev nodes

ErrorPlot(E,K,99)

% =========================================================================

function [p,err] = DrawInterp(f,a,b,n,type,sing,sub)

% parameter --------------------------------------------------------

N = 2000;                       % number of points for plotting/error
n_max = 125;                    % maximum degree to be plotted
y_min = -0.5; y_max = 1.5;      % y range
x_str = -0.2; y_str = 1.125;    % position of text

% ------------------------------------------------------------------

x = linspace(a,b,N);

switch type
    case 'equi'                             % equidistant nodes
        knots = linspace(a,b,n+1);
        title_str='equidistant Nodes';
    case 'cheb'                             % Chebyshev nodes
        knots = (a+b)/2+(b-a)/2*cos(linspace(0,pi,n+1));
        title_str='Chebyshev Nodes';
end
lambda = weights(knots); fval = f(knots);
p = barycentric1(x,knots,fval,lambda);

if n <= n_max
    subplot(1,2,sub);
    h0 = plot(x,p,'b-',knots,fval,'r.'); hold on;
    h1 = plot([-sing,-sing],[y_min,y_max],'g-',...
        [sing,sing],[y_min,y_max],'g-'); hold off;
    set(h0,'LineWidth',2,'MarkerSize',16);
    set(h1,'LineWidth',5);
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

function ErrorPlot(E,K,fig)

m = length(E);
h=figure(fig); set(h,'Name','error plot')
h=semilogy(E,'b-'); hold on
set(h,'LineWidth',2,'MarkerSize',16);
p = polyfit(1:20,log(E(1:20)),1);
C= exp(p(2));
h=semilogy(1:m,C./K.^(1:m),'r-');
set(h,'LineWidth',2,'MarkerSize',16);
axis([0 m, 10^-15, 10^0])
title('error of interpolation in in Chebyshev nodes')
xlabel('polynomial degree');
ylabel('absolute error');


% =========================================================================

function fig = initdraw

close all
screensize=get(0,'ScreenSize');
fig=figure('Outerposition',screensize,'Units',get(0,'Units'));
figure(fig);
set(fig,'DoubleBuffer','on','Name','Runge-Phänomen');

