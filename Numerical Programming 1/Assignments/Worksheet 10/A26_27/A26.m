function A26

n=10;
A=randn(1000);
[mu,omega]=Rayleigh(A,rand(size(A,2),1),n);


semilogy(1:n,max(omega,eps),'.');
xlabel('iteration number');
ylabel('backward error');
display(mu);

end