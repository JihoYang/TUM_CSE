m = 200; n = 10;
A = randn(m,n);
spy(A), shg
B = A'*A;
spy(B), shg
spy(A*A'), shg
R = chol(B);
spy(R), shg
Q = A*inv(R);
norm(A-Q*R)

ans =

      2.21932948609256e-15

norm(Q'*Q-eye(n))

ans =

      4.27529620311231e-16

spy(Q), shg
m = 2000; n = 1000;
A = randn(m,n);
B = A'*A;
R = chol(B);
Q = A*inv(R);
norm(A-Q*R)

ans =

      4.94446206119011e-14

norm(Q'*Q-eye(n))

ans =

      3.46506746158196e-15

m = 20; n = 10;
A = randn(m,n);
A = triu(A);
B = A'*A;
R = chol(B);
Q = A*inv(R);
norm(Q'*Q-eye(n))

ans =

      7.71464907185155e-09

condest(B)

ans =

          19915326228.4689

E = Q'*Q;
E(10,10)

ans =

          1.00000000770107

[Q,R] = qr(A,0);
spy(Q), shg
norm(Q'*Q-eye(n))

ans =

     0

m = 2000; n = 1000;
A = randn(m,n);
A = triu(A);
R = chol(B);
B = A'*A;
[Q,R] = qr(A,0);
spy(Q), shg
A = randn(m,n);
[Q,R] = qr(A,0);
spy(Q), shg
[Q,R] = qr(A);
spy(Q), shg
size(Q)

ans =

        2000        2000

[Q,R] = qr(A,0);
size(Q)

ans =

        2000        1000

n = 100;
A = randn(n,n);
b = randn(n,1);
xx = A\b;
[Q,R] = qr(A);
x = R\(Q'*b);
norm(x-xx,inf)

ans =

      5.82645043323282e-13

r = b - A*x;
rr = b - A*xx;
omega = norm(r,inf)/(norm(A,inf)*norm(x,inf))

omega =

      3.76180258190835e-17

omega_xx = norm(rr,inf)/(norm(A,inf)*norm(xx,inf))

omega_xx =

      6.36833935291599e-17

n = 4000;
A = randn(n,n);
tic, [L,U,p] = lu(A,'vector'); toc
Elapsed time is 0.701033 seconds.
tic, [Q,R] = qr(A); toc
Elapsed time is 3.249357 seconds.
