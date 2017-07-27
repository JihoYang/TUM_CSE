format short g

N = 1e6;
X = rand(N,1);
f = @(x) sin(pi*x);

[mean(f(X));2/pi]

ans =

      0.63697
      0.63662

N = 1e8;
X = rand(N,1);
tic, X = rand(N,1); toc
Elapsed time is 1.274951 seconds.
tic, Q = mean(f(X)), toc

Q =

      0.63663

Elapsed time is 0.798248 seconds.
[Q;2/pi]

ans =

      0.63663
      0.63662

tic, X = randn(N,1); toc
f = @(x) x.^2;
mean(f(X))

ans =

       1.0001

% example with normal distribution on R^d

d = 10;
f = @(x) sum(x.^2)./(1+sum(x.^2));

N = 1e5;
X = randn(d,N);
fX = bsxfun(@(x,y) f(x), X, 0);
size(fX)

ans =

           1      100000

Q = mean(fX)

Q =

      0.89182

I = 0.8920274601139534

I =

      0.89203

d = 1000;
tic, X = randn(d,N); toc
Elapsed time is 1.712275 seconds.
fX = bsxfun(@(x,y) f(x), X, 0);
size(fX)

ans =

           1      100000

Q = mean(fX)

Q =

        0.999

I = 0.9989990010090231

I =

        0.999

err = abs(Q-I)

err =

    5.791e-08

format long g
[Q;I]

ans =

         0.998999058919046
         0.998999001009023

d = 1000;
N = 1e5;
X = randn(d,N);
X = rand(d,N);
f = @(x) sum(x.^2);
fX = bsxfun(@(x,y) f(x), X, 0);
size(fX)

ans =

           1      100000

Q = mean(fX)

Q =

          333.328565865804

