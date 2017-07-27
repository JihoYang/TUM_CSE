function A15

f = @(t) t.^2 .* (t-1).^2 .*(t-1/200).^2 .* (t-1/100).^2 .* (t-101/200).^2;

TOLs=10.^(-2:-1:-8);

for tol=TOLs
    result=quad24(f,0,1,tol);
    fprintf('tol=%4.1e, result=%e\n',tol,result);
end



end