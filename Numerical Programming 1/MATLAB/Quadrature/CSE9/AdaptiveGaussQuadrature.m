function [val,err,nfcn] = AdaptiveGaussQuadrature(f,interval,TOL,m)

% function [val,err,nfcn] = AdaptiveGaussQuadrature(f,interval,TOL,m)
%
% Lecture "Numerische Mathematik II"
% Folkmar Bornemann TUM, SoSe 2002
%
% idea following E. Hairer, Introduction à l'Analyse Numérique, 
%     Université de Genève, 2001, S. 17-18
%
% adaptive Gaussian quadrature with step-size control
%					                   
% input:        f	        function
%               interval    interval of form [a b]
%               TOL 	    relative user tolerance
%		        m	        p = 8m+6 order of method
%
% output:       val	        computed value of integral
%		        err	        estimate of relative error
%		        nfcn	    number of f-evaluations
%
% example:
%		f = @(x) 1./(1e-4+x.^2);
%		[val,err,nfcn] = AdaptiveGaussQuadrature(f,[-1 1],1e-11,2)

global nfcn w dw1 dw2 c

p = embedded(m); nfcn = 0; 
a = interval(1); b= interval(2);

x = a; h = (b-a); result = []; val = 0; err = 0;
val_scl = abs(quadrature(f,a,h));
TOL = TOL*val_scl;

while x(end) < b
  while 1
    [val_loc,err_loc] = quadrature(f,x(end),h(end));
    h_new = min(2,0.75*(TOL/err_loc)^(1/(p+1)))*h(end);
    if err_loc <= TOL, break; end
    h(end) = h_new;
  end
  x(end+1) = x(end) + h(end);
  h(end+1) = min(h_new,b-x(end));
  val = val + val_loc;
  err = err + err_loc;
end
graphics(f,x,h);

err = err/abs(val);

return

%-----------------------------------------------------------------------------------

function [val,err] = quadrature(f,x,h)

global nfcn w dw1 dw2 c;

fval = f(x+h*c); nfcn = nfcn + length(c);
err1 = h*dw1*fval; 
err2 = h*dw2*fval;
val = h*w*fval;
err = abs(err1)*(err1/err2)^2;

return;

%-----------------------------------------------------------------------------------

function p = embedded(m); % not not look at this, it could be a table

global w dw1 dw2 c;

s = 4*m+3; s1 = 4*m+2; s2 = 2*m;
w1 = zeros(1,s); w2 = zeros(1,s);
e1 = zeros(s,1); e1(1) = 1;

ind1 = 1:s; ind2 = ind1; 
ind2(1:2:s) = []; ind2(m+1) = [];
ind1(2*m+2) = []; 

[w,c,P] = gauss_legendre(s);
k1 = length(ind1); w1(ind1) = P(1:k1,ind1)\e1(1:k1);
k2 = length(ind2); w2(ind2) = P(1:k2,ind2)\e1(1:k2);

dw1 = w-w1; dw2 = w-w2;
p = 2*s;

return

%-----------------------------------------------------------------------------------

function [w,c,P] = gauss_legendre(s) % not not look at this, it could be a table

k = 1:s-1; beta = k./sqrt((2*k-1).*(2*k+1));
T = diag(beta,-1) + diag(beta,1);
[Q,L] = eig(T);
c = (diag(L)+1)/2; 
w = Q(1,:).^2; 
P = Q*diag(1./Q(1,:));

return;

%-----------------------------------------------------------------------------------

function graphics(f,x,h)

global nfcn w; warning off;

s = length(w);
a = x(1); b = x(end);
x_plot = linspace(a,b,1000);
f_plot = f(x_plot); 
mn = min(f_plot); mx = max(f_plot);
mf1 = mn-(mx-mn)/20;
mf2 = mn-2*(mx-mn)/20;
mf3 = mx+2*(mx-mn)/20;

figure('Name','adaptive Gauss quadrature with step size control');
nin = length(x)-1; ovh = 100*(nfcn - s*nin)/s/nin;
str_title = sprintf('order %d ---- statistics: %d intervals, %2.0f%% overhead for adaptivity',2*s,nin,ovh);
plot(x_plot,f_plot,'g-','LineWidth',2,'LineSmoothing','On');
hold on;
plot(x,mf1,'r+',[a b],[mf1 mf1],'k-',[0 0],[mf2 mf3],'k:',...
  x,f(x),'r.',[a b],[0,0],'k:');
axis([a b mf2 mf3]);
title(str_title);
xlabel('x');
ylabel('f');

warning on;

return;