function val = trapeziodal(f,x) % composite trapezoidal rule with nodes x

y = f(x);
val = dot(diff(x),(y(1:end-1)+y(2:end))/2);
