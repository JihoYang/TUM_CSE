function A2
% Assignment 2

f=@(x) (1-cos(x))./x;

j=4:9;x = 10.^(-j);
erg=[x; f(x); sol(x)];

fprintf('\n\n');
fprintf('    x                   | (1-cos(x))/x          |');
fprintf('  (2*sin(x/2)^2)/x     |\n');
fprintf('=================================================');
fprintf('========================\n');
fprintf(' %1.15e  | %1.15e | %1.15e |\n',erg);

xx=linspace(1e-9,1e-7,5000);
plot(xx,f(xx),xx,sol(xx));
legend({'(1-cos(x))/x','(2*sin(x/2)^2)/x'},'Location','NorthWest');

end

function erg=sol(x)
erg = (2*sin(x/2).^2)./x;
erg(x==0)=0;
end

