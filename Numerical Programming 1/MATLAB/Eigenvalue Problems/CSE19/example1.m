lambda=1:25;
S=rand(numel(lambda));
A=S\(diag(lambda)*S);

B=qr_basic(A,150);
