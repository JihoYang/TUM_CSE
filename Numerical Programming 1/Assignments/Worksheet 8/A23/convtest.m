function convtest
% some test examples

rng(0815); n=1e3; A=sprandsym(n,0.1,1e-2,2); 
figure(1);convplot(A,ones(n,1),250);

rng(0815); n=3e3; A=sprandsym(n,0.1,1e-2,2); 
figure(2);convplot(A,ones(n,1),250);

rng(0817); n=1e3; A=sprandsym(n,0.1,1e-3,2);
figure(3);convplot(A,ones(n,1),3000);

m=200; A=laplace2DuniformA(m); n=m*m;
figure(4);convplot(A,ones(n,1),400);

n=1e4; A=spdiags([4*ones(n,1),-2*ones(n,1),ones(n,1)],[0,-1,1],n,n);
figure(5);convplot(A,A*ones(n,1),100);

end

function A=laplace2DuniformA(n)
I=speye(n,n); one = ones(n,1);
D=spdiags( [-4*one, one, one],[0,-1,1], n,n);
C=spdiags( [one,one],[-1,1], n,n);
A=(kron(I,D) + kron(C,I))*(n+1)*(n+1);
end
