function pts=A6(a,b,n)
% generate n Chebyshev points in interval (a,b).
assert( a<=b );
nm1 = n-1;
x = cos( pi/nm1*(nm1:-1:0) );            % Chebyshev points in (-1,1)
if a>0 || b<0
    len_half = (b-a)/2;
    pts = x*len_half + (a+len_half);
else
    m = (a+b)/2;
    pts = x*(b-m)+m;
end
end