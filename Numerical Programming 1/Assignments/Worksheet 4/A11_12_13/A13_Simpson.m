function result = A13_Simpson(a,b,f,n)

assert(n>1);
x=linspace(a,b,2*n-1)';   % simpson nodes; x_j, (x_j+x_(j+1))/2, x_(j+1)

w=repmat([2,4],1,n);                 % weights: [1,4,2,4,2,4,...,4,1]*h/6
w(1)=1;w(2*n-1)=1; w(2*n:end)=[];
w=w*(b-a)/(6*(n-1));

result = w*f(x);

end