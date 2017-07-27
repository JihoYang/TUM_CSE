e = 1e-7;
A = [1 1; e 0; 0 e]; b = [2;e;e]; % LŠuchli example
x = (A'*A)\(A'*b) % solution of normal equation

x =

          1.01123595505618
          0.98876404494382

e = 1e-8;
A = [1 1; e 0; 0 e]; b = [2;e;e]; % LŠuchli example
x = (A'*A)\(A'*b) % solution of normal equation
[Warning: Matrix is singular to working precision.] 

x =

   NaN
   NaN

e = 1e-7;
A = [1 1; e 0; 0 e]; b = [2;e;e]; % LŠuchli example
[Q,R] = qr(A);
x = R\(Q'*b)

x =

                         1
                         1

e = 1e-8;
A = [1 1; e 0; 0 e]; b = [2;e;e]; % LŠuchli example
[Q,R] = qr(A);
x = R\(Q'*b)

x =

                         1
                         1

m = 8000; n = 2000;
A = randn(m,n);
b = randn(m,1);
tic, [Q,R] = qr(A,0); toc
Elapsed time is 1.980790 seconds.
tic, x = R\(Q'*b); toc
Elapsed time is 0.018426 seconds.
tic, xx = A\b; toc
Elapsed time is 7.452752 seconds.
norm(x-xx,2)

ans =

      2.29020779409654e-15

rho = norm(b-A*x,2)

rho =

          77.9249591514244


tic, R1 = triu(qr([A,b])); toc
Elapsed time is 1.034934 seconds.
tic, x1 = R1(1:n,1:n)\R1(1:n,n+1); toc
Elapsed time is 0.023266 seconds.
norm(x-x1,2)

ans =

      1.20158779796254e-15

[rho abs(R1(n+1,n+1))]

ans =

          77.9249591514244          77.9249591514242

