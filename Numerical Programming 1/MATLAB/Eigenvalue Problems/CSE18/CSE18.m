
A = diag(1:22)

chi = charpoly(A);


lambda = roots(chi);
lambda(7:8)

ans =

           15.4541352582099 +      1.17490829850328i
           15.4541352582099 -      1.17490829850328i

plot(lambda,'b*'), shg

lambda = 500*0.8.^(1:1000);
lambda(1)

ans =

   400

lambda(1000)

ans =

      6.15115961080593e-95

[Q,~] = qr(randn(1000,1000));
A = Q'*diag(lambda)*Q;
format short g; A(1:5,1:5)

ans =

       4.6886      -0.1181     -0.11663       1.7336      -1.9887
      -0.1181       1.0597      0.24493    -0.044642     -0.13186
     -0.11663      0.24493       1.3693     -0.58548     0.058988
       1.7336    -0.044642     -0.58548       2.3541      -0.9805
      -1.9887     -0.13186     0.058988      -0.9805       2.5931

tic, [mu,v,omega,k] = PowerMethod(A,rand(1000,1),1e-12); toc
Elapsed time is 0.061082 seconds.
omega

omega =

   8.1225e-13

mu

mu =

   400

lambda = -500*0.8.^(1:1000);
A = Q'*diag(lambda)*Q;
tic, [mu,v,omega,k] = PowerMethod(A,rand(1000,1),1e-12); toc
Elapsed time is 0.057139 seconds.
mu

mu =

         -400

k

k =

   114

log(1e-12)/log(0.8)

ans =

       123.83

0.8^124

ans =

   9.6196e-13

lambda = -500*0.95.^(1:1000);
A = Q'*diag(lambda)*Q;
tic, [mu,v,omega,k] = PowerMethod(A,rand(1000,1),1e-12); toc
Elapsed time is 0.212593 seconds.
k

k =

   491

log(1e-12)/log(0.95)

ans =

       538.69

abs(lambda(2))/abs(lambda(1))

ans =

         0.95

lambda = -500*0.95.^(1:1000);
lambda = [lambda 500*0.95.^(1:1000)];
[Q,~] = qr(randn(2000,2000));
A = Q'*diag(lambda)*Q;
tic, [mu,v,omega,k] = PowerMethod(A,rand(2000,1),1e-12); toc
{Operation terminated by user} 

lambda = 500*0.95.^(1:1000);

lambda(85)

ans =

       6.3896

[Q,~] = qr(randn(1000,1000));
A = Q'*diag(lambda)*Q;
tic, [mu,v,omega,k] = InversePowerMethod(A,6,rand(1000,1),1e-12); toc
Elapsed time is 0.072571 seconds.
mu

mu =

       6.0702

lambda(86)

ans =

       6.0702

lambda(87)

ans =

       5.7667

theta = abs((6.0702-6)/(5.7667-6))

theta =

       0.3009

log(1e-12)/log(theta)

ans =

       23.007

k

k =

    15

lambda = -500*0.95.^(1:1000);
lambda = [lambda 500*0.95.^(1:1000)];
[Q,~] = qr(randn(2000,2000));
A = Q'*diag(lambda)*Q;
tic, [mu,v,omega,k] = InversePowerMethod(A,lambda(1),rand(2000,1),1e-12); toc
Elapsed time is 0.136223 seconds.
mu

mu =

         -475

lambda(1)

ans =

  -475

k

k =

     1

tic, [mu,v,omega,k] = InversePowerMethod(A,lambda(1)+eps,rand(2000,1),1e-12); toc
Elapsed time is 0.135081 seconds.
k

k =

     1