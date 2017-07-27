
n=150; prange=[1,2,3,4,inf];
figure(1); cm=colormap('lines');
for k=1:numel(prange)
    p=prange(k);
    c=A17_contour(p,n);hold on;
    set(c,'LineColor',cm(k,:));
end
hold off;axis equal;
legend(arrayfun( @(p) ['p=',num2str(p)], prange, 'UniformOutput',0), ...
    'Location','Best');
