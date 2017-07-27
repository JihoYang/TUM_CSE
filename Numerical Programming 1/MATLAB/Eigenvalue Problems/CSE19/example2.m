lambda=1:25;
S=rand(numel(lambda));
A=S\(diag(lambda)*S);

[B,iter]=qr_iteration(A);
display(iter);
