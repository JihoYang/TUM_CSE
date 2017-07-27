function [mu,omega]=Rayleigh (A,u,n)
% n steps of the Rayleigh iteration 

omega=NaN*ones(1,n); mu=omega;
normA=norm(A,'fro');
I=eye(size(A));

for k=1:n
    u=u/norm(u,2);
    Au=A*u;
    mu(k) = u'*Au;
    omega(k)=norm(mu(k)*u - Au,2)/normA;
    u = (A-mu(k)*I)\u;
end

end