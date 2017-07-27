function pp = barycentric2(xx,x,f,lambda)

% function pp = barycentric2(xx,x,f,lambda)
%
% Auswertung der 2. baryzentrischen Interpolationsformel
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

numer = zeros(size(xx)); 
denom = zeros(size(xx)); 
exact = zeros(size(xx));

for j = 1:length(x)
    xdiff = xx-x(j);
    temp  = lambda(j)./xdiff;
    numer = numer + temp*f(j);
    denom = denom + temp;
    exact(xdiff==0) = j;
end

pp = numer./denom;
jj = find(exact); pp(jj) = f(exact(jj));
