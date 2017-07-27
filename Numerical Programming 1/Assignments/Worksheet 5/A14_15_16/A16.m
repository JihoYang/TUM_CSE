function A16

a=0; b=1; exact = -4/9; Iscal=1; f = @sqrtlog;

% the algorithms we are testing
%       [val,err,nfcn] = algorithm(f,a,b,tol,Iscal)
algorithms=struct(...
    'name', {'quad24','A14','adapt. Gauss (m=0)',...
             'adapt. Gauss (m=1)','adapt. Gauss (m=2)' },...
    'func', {
    @(f,a,b,tol,Iscal) quad24(f,a,b,tol*Iscal),...
    @(f,a,b,tol,Iscal) A14(f,a,b,tol,Iscal),...
    @(f,a,b,tol,Iscal) AdaptiveGaussQuadrature(f,[a,b],tol*Iscal,0,1),...
    @(f,a,b,tol,Iscal) AdaptiveGaussQuadrature(f,[a,b],tol*Iscal,1,1),...
    @(f,a,b,tol,Iscal) AdaptiveGaussQuadrature(f,[a,b],tol*Iscal,2,1) } );


TOLs=kron(10.^(-2:-1:-15),[1,0.5]);  % the tolerances
TOLs(end)=[];

no_algorithms=numel(algorithms);    %  number of algorithms
no_TOLs=numel(TOLs);                %  number of tolerances


%% Calculate
val=NaN*ones(no_TOLs,no_algorithms); err=val; nfcn=val;
for alg=1:no_algorithms
    fprintf('  Algorithmn %s TOLs:',algorithms(alg).name);
    for tol=1:no_TOLs
        fprintf(' %g',TOLs(tol));
        [val(tol,alg),err(tol,alg),nfcn(tol,alg)]=...
                algorithms(alg).func(f,a,b,TOLs(tol),Iscal);
            err(tol,alg) = err(tol,alg)/abs(val(tol,alg)); % rel. error  
    end    
    fprintf('\n');
end
    

%% Output:
figure; subplot(2,1,1);
loglog(nfcn,max(abs( (val-exact)./exact),eps));
xlabel('#f-Evaluations');ylabel('rel. error');
legend(algorithms(:).name);
grid on;

subplot(2,1,2);
h = max(abs( (val-exact)),eps);loglog( h, err );hold on;  
h = [ min(h(:)),max(h(:)) ];loglog(h,h,'k--');hold off;
xlabel('error');ylabel('estimated error');
grid on;

end

function y=sqrtlog(x)
y=sqrt(x).*log(x);
y(x==0)=0;
end
