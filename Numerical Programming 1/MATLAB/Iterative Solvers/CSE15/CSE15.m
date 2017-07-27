
	Academic License

n = 1000;
Generate
whos
  Name         Size               Bytes  Class     Attributes

  A         1000x1000            135488  double    sparse    
  b         1000x1                 8000  double              
  n            1x1                    8  double              
  rho          1x1                    8  double              

spy(A), shg
B = full(A);
whos
  Name         Size                Bytes  Class     Attributes

  A         1000x1000             135488  double    sparse    
  B         1000x1000            8000000  double              
  b         1000x1                  8000  double              
  n            1x1                     8  double              
  rho          1x1                     8  double              

nnz(A)/n^2

ans =

                  0.010957

tic, [L,U,p] = lu(A); toc
Elapsed time is 0.363991 seconds.
spy(L); shg
figure(2); spy(A); shg
nnz(L)/n^2/2

ans =

                  0.163568

n = 10000;
Generate
whos
  Name          Size                 Bytes  Class     Attributes

  A         10000x10000            1359380  double    sparse    
  B          1000x1000             8000000  double              
  L          1000x1000             3964712  double    sparse    
  U          1000x1000             3964712  double    sparse    
  ans           1x1                      8  double              
  b         10000x1                  80000  double              
  n             1x1                      8  double              
  p          1000x1000               16004  double    sparse    
  rho           1x1                      8  double              

1359380/(8*n^2)

ans =

               0.001699225

nnz(A)/n^2

ans =

                0.00109948

tic, [L,U,p] = lu(A); toc
Elapsed time is 349.886596 seconds.

350/60

ans =

          5.83333333333333

x = randn(n,1);
tic, y = A*x; toc
Elapsed time is 0.002091 seconds.
B = full(A);
tic, y = B*x; toc
Elapsed time is 0.104318 seconds.
0.1/0.002

ans =

    50

1/(nnz(A)/n^2)

ans =

           909.52086440863

b = randn(n,1);
tic, z = L\b(p); toc
{Out of memory. Type HELP MEMORY for your options.
} 

StationaryIteration
Elapsed time is 0.082258 seconds.

ans =

                        88       8.4000639350195e-11

n = 100000;
Generate
StationaryIteration
Elapsed time is 0.942058 seconds.

ans =

                        89       9.1896724245694e-11

% now change to Gauss-Seidel!

StationaryIteration
Elapsed time is 0.573870 seconds.

ans =

                        49      7.37351563907342e-11



