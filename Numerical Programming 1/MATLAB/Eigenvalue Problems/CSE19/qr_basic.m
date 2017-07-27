function A=qr_basic(A,n)

opt=struct('levels',10.^(-15:2:2));
res=showMatrixAbs('init',opt);showMatrixAbs('update',res,A);

for k=1:n
    figure(res.fig);title(res.ax1,['Iteration ',num2str(k)]);drawnow;
    [Q,R]=qr(A);                     % QR-decomposition
    A = R*Q;                         % next iterate
    showMatrixAbs('update',res,A);   % update figure
    display(diag(A).');               % print out diagonal
end

end