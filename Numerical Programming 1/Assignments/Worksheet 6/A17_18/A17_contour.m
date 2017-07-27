function c=A17_contour(p,n)
% Visualize set { ||x||_p = 1 }
% Input:
%     p         p for Norm (typical: p>=1)
%     n         number of nodes for (visualization) grid

[X,Y]=meshgrid(  2*(-1:n+1)/n -1 );
if p<inf
    Z=abs(X).^p + abs(Y).^p;
else
    Z=max(abs(X),abs(Y));
end

[~,c]=contour(X,Y,Z,[1,1]);

end