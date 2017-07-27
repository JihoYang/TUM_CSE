function err = A13_err(a,b,f,nrange,exact)

err=zeros(1,numel(nrange));
for k=1:numel(nrange)
    value = A13_Simpson(a,b,f,nrange(k));
    err(k) = max(abs(value-exact),eps);
end

end