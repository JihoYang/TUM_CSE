function [x,fail] = Newton(f,Df,x,tol)

y = f(x); J = Df(x);
while true		
    [L,U] = lu(J);                         
    dx = -U\(L\y); dx_norm = norm(dx);     
    if isnan(dx_norm), fail = true; break; end
    x = x + dx;
    if dx_norm <= tol, fail = false; break; end   
    y = f(x); J = Df(x);             
    dx_bar = norm(-U\(L\y));
    if dx_bar > dx_norm, fail = true; break; end
end