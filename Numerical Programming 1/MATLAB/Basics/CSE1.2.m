%-- 14.10.15, 15:00 --%
x = 1;
e = 0;
while x<inf, x=2*x; e=e+1; disp([x,e]); end

x = single(1);
e = 0;
while x<inf, x=2*x; e=e+1; disp([x,e]); end

realmax
realmax('single')

x = 1;
e = 0;
while x>0, x=x/2; e=e-1; disp([x,e]); end

x = single(1);
e = 0;
while x>0, x=x/2; e=e-1; disp([x,e]); end

realmin
realmin('single')

(1e-16 + 1) - 1
1e-16 + (1 - 1)
(eps + 1) - 1
(eps/2 + 1) - 1

x = single(0.1)
num2hex(x)
(0.1 -x)/0.1