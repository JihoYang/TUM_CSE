function A10b(a,b,f,fP,nrange)
% convergence plots
%
% examples:
% A10b(0,2*pi,@sin,@cos,round(logspace(1,3,15)));
%
% f = @(x) (x.^2 + 10./(sin(x)+1.2))/45
% fP = @(x)  (x-5*cos(x)./((6/5+sin(x)).^2))*2/45
% A10b(-2,8,f,fP,round(logspace(1,3,30)));

nrange=sort(nrange);

%% Show figure with graphfs for smallest nrange(1)
x=linspace(a,b,nrange(1)); y=f(x);
xx=linspace(a,b,1e3);
[yy_spline,yy_hermite]=A10a(x,y,fP(x),xx);

figure(1);
plot(x,y,'r*',xx,f(xx),'k',xx,yy_spline,'b',xx,yy_hermite,'g');
legend({'nodes: x_j','f','spline','hermite'});
title([func2str(f),' with ',num2str(nrange(1)),' nodes ']);

%% Error plots
err = zeros(2,numel(nrange)); % err_spline and err_hermite

for k=1:numel(nrange)
    x=linspace(a,b,nrange(k));
    xx=linspace(a,b,max(500,10*nrange(k)));
    yy=f(xx);
    [yy_spline,yy_hermite]=A10a(x,f(x),fP(x),xx);
    err(1,k) = max(eps,norm(yy-yy_spline,inf));
    err(2,k) = max(eps,norm(yy-yy_hermite,inf));
end

figure(2);
loglog(nrange,err,'.','MarkerSize',18);xlabel('n');ylabel('error');
title(['convergence plots for ',func2str(f)]);
grid on;legend({'spline','hermite'});

end