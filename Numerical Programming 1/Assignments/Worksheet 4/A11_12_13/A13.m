function A13

figure(1);
f = @(x) cosh(sin(x)); exact = 2.283194520910577988;
%nrange = round(logspace(1,4,40));
nrange = 2:17;
err = A13_err(-1,1,f,nrange,exact);
semilogy(nrange,err,'.','MarkerSize',20);grid on;
xlabel('n'); ylabel('error');
title(['convergence plot for ',func2str(f)]);

figure(2);
f = @(x) abs(sin(20*x)); exact = (6+sin(10)^2)/5;
nrange = round(logspace(1,6,40));
err = A13_err(-1,1,f,nrange,exact);
loglog(nrange,err,'.','MarkerSize',20);grid on;
xlabel('n'); ylabel('error');
title(['convergence plot for ',func2str(f)]);

shg;

end

