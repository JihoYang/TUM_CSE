warning off symbolic:sym:int:notFound

syms x h;
f=symfun(sym('f(x)'),x);

q1=@(h) h*( f(0)+3*f(h/3)+3*f(2*h/3)+f(h) )/8;
E1=@(h) int(f(x),x,0,h)-q1(h);
disp('Error term (a)');
pretty(simplify(taylor(E1(h),h,'Order',6)))

alpha=0*h+sqrt(15)/10;
q2=@(h) (h/18)*( 5*f((-alpha+1/2)*h) + 8*f(h/2) + 5*f((alpha+1/2)*h));
E2=@(h) int(f(x),x,0,h)-q2(h);
disp('Error term (b)');
pretty(simplify(taylor(E2(h),h,'Order',8)))

