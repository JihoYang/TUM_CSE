function A18
f = @cos; I_exact = sqrt(2*pi/exp(1));         % test-function & exact
N=1e3; krange=unique(round(logspace(0,6,40))); % values of N and k
val = zeros(1,numel(krange));                  % memory to save result

for i=1:numel(krange)
    k=krange(i);
    fprintf('k=%8i, N=%i\n',k,N);drawnow;
    if i>1 % use same trick als in mc_composite and simulate only ...
        kOld=krange(i-1); d=k-kOld; % ... d*N new experiments
        val(i) = (val(i-1)*kOld + mc_composite(f,N,d)*d)/k;
    else
        val(i) = mc_composite(f,N,k);
    end
end
err = max(eps,abs((val-I_exact)./I_exact));
loglog(N*krange,err,'.');xlabel('number of nodes');ylabel('rel. error');
end

function result = mc_composite(f,N,k)
% A18 (c); memory: O(N)
result = 0;
for j=1:k                       % k calls to mc
    result = result + mc(f,N);  % sum all values
end
result = result/k;              % calculate mean
end

function result = mc(f,N)
% A18 (a): normalizing factor is 1/sqrt(2*pi); memory: O(N)
result = sqrt(2*pi)*mean(f(randn(N,1)));
end