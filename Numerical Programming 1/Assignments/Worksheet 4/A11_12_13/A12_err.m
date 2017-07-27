function err = A12_err(a,b,f,nrange,exact)

err=zeros(1,numel(nrange));
for k=1:numel(nrange)
    [x,w] = legpts(nrange(k),[a,b]);
    err(k) = max(abs(w*f(x)-exact),eps);
end

end