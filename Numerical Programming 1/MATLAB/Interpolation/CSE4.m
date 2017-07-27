
	Academic License

n=10; omega = exp(-2*pi*1i/n);
M = zeros(n,n);
for j=0:n-1, for k=0:n-1, M(j+1,k+1) = omega^(-j*k); end, end
format short g
M(1:3,1:3)

ans =

            1 +          0i            1 +          0i            1 +          0i
            1 +          0i      0.80902 +    0.58779i      0.30902 +    0.95106i
            1 +          0i      0.30902 +    0.95106i     -0.80902 +    0.58779i

U = M/sqrt(n);
norm(U*U' - eye(n))

ans =

   6.2118e-15

n = 2000; omega = exp(-2*pi*1i/n);
M = zeros(n,n);
for j=0:n-1, for k=0:n-1, M(j+1,k+1) = omega^(-j*k); end, end
y = randn(n,1);
tic, c = (1/n)*M'*y; toc
Elapsed time is 0.072053 seconds.


n = 1e7;
time = (n/2000)^2*0.072

time =

      1.8e+06

time/3600/24

ans =

       20.833

y = randn(n,1);
tic, c = (1/n)*ifft(y); toc
Elapsed time is 0.429444 seconds.

n=1e7

n =

    10000000

8*n^2

ans =

        8e+14

5*n*log2(n)

ans =

   1.1627e+09

help fft

doc fft
