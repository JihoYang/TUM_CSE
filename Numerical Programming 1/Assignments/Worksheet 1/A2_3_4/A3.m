function A3
% Assignment 3

f=@(x) sqrt(x+1)-sqrt(x);
solf = @(x) 1./(sqrt(x+1)+sqrt(x));

g=@(x) asin(2*(1-cos(x))./(x.^2));
solg = @(x) asin( (2*sin(x/2)./x).^2 );

x = 10.^(48:53);
erg=[x; f(x); solf(x)];

fprintf('\n\n');
fprintf('    x                   | sqrt(x+1)-sqrt(x)     |');
fprintf('1/(sqrt(x+1)+sqrt(x))  |\n');
fprintf('=================================================');
fprintf('========================\n');
fprintf(' %1.15e  | %1.15e | %1.15e |\n',erg);

x = 10.^(-(4:10));
erg=[x; g(x); solg(x)];

fprintf('\n\n');
fprintf('    x                   |asin(2*(1-cos(x))/(x^2))|');
fprintf('asin((2*sin(x/2)./x)^2)|\n');
fprintf('=================================================');
fprintf('========================\n');
fprintf(' %1.15e  | %1.15e | %1.15e |\n',erg);

  
end

% function erg=solg2(x)
% w=cos(x);
% erg= asin( 2*(1-w)./(acos(w)).^2 );
% erg(w==1)=pi/2;
% end


