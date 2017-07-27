function [A,iter]=qr_iteration(A)

opt=struct('levels',10.^(-15:2:2));
res=showMatrixAbs('init',opt);showMatrixAbs('update',res,A);
m=size(A,1); % working on rows/columns 1:m
I=eye(m); normA=norm(A,'fro');

A=hess(A);                    % preprocessing: hessenberg form
iter=zeros(1,m);

while (m>1)    
    figure(res.fig);          % drawing
    title(res.ax1,['Iteration ',num2str(iter(m)),' for m=',num2str(m)]);
    drawnow;
    
    mu = wilkinson_shift(A(m-1:m,m-1:m));    % Shift
    [Q,R]=qr(A(1:m,1:m)-mu*I);               % QR-decomposition with shift
    A(1:m,1:m) = R*Q+mu*I;                   % next iterate in 1:m block
    iter(m)=iter(m)+1;                       % increate iteration counter
    omega = norm(A(m,1:m-1))/normA;          % backward error
    if omega<1e-15                           % backward stable?
        fprintf(['found backward stable eigenvalue %f; ',...
            'now working on first %d rows/columns\n'],A(m,m),m-1);
        A(m,1:m-1)=0;                        % numerical deflation
        m=m-1; I=eye(m);                     % work only on smaller block
    end
    
    showMatrixAbs('update',res,A);           % update figure
end

end


%% choice of shift

function mu=wilkinson_shift(B)
cand = eig(B);
if abs(cand(1)-B(end,end))<abs(cand(2)-B(end,end))
    mu=cand(1);
else
    mu=cand(2);
end
end