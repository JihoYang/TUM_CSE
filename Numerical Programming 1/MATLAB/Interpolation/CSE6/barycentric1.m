function pp = barycentric1(xx,x,f,lambda)

% function pp = barycentric1(xx,x,f,lambda)
%
% Auswertung der 1. baryzentrischen Interpolationsformel
% (c) F. Bornemann 24/01/04
%
% Input:
%           xx      Argument des Interpolationspolynoms (Vektor)
%           x       Stützstellen (Vektor)
%           f       Stützwerte   (Vektor)
%           lambda  Gewichte     (Vektor)
%
% Output:
%           pp      Wert des Interpolationspolynoms

sum   = zeros(size(xx)); 
omega = ones(size(xx)); 
exact = zeros(size(xx));

for j = 1:length(x)
    xdiff = xx-x(j);
    sum   = sum + lambda(j)*f(j)./xdiff;
    omega = omega.*xdiff;
    exact(xdiff==0) = j;
end

pp = omega.*sum;
jj = find(exact); pp(jj) = f(exact(jj));
