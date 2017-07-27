function lambda = weights(x)

% function lambda = weights(x)
%
% Rekursive Berechnung der Gewichte für die Lagrange-Interpolation
% (c) F. Bornemann 24/01/04
%
% Input:
%           x       Stützstellen
%
% Output:
%           lambda  Gewichte  

lambda = ones(size(x));

for n=1:length(x)
    for j=1:n-1
        d = x(j)-x(n);
        lambda(j) =  lambda(j)/d;
        lambda(n) = -lambda(n)/d;
    end
end
