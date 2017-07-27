function [E,C,residual]=A21_solve(filename)
% transform, setup linear least square problem and solve
% Input:
%   filename        file with measurements

M=load(filename);m=size(M,1);n=2;  % T_j, K_j, ΔK_j
D=spdiags(M(:,2)./M(:,3),0,m,m);   % K_j/ΔK_j

A=D*[-1./M(:,1),ones(m,1)];        % 1. column -1/T_i, 2. column: 1
b=D*log(M(:,2));

R1=triu(qr([A,b]));                % use "Q-free" method to solve ...
x=R1(1:n,1:n)\R1(1:n,n+1);         % ... least square problem
E=x(1); C=exp(x(2));               % calculate E and C from x
residual=R1(n+1,n+1);
end