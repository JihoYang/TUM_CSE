function X=A27(A,B,C)
n=length(A);
[U,R]=schur(A,'complex');            % O(n^3)
[V,S]=schur(B,'complex');            % O(n^3)
E=U'*C*V;                            % O(n^3)
Y=zeros(n); I=eye(n);
for j=1:n                            %    n times  ...
  Y(:,j)=(R-S(j,j)*I)\(E(:,j) + Y(:,1:j-1)*S(1:j-1,j));   % ... O(n^2)
end
X=U*Y*V';                            % O(n^3)
