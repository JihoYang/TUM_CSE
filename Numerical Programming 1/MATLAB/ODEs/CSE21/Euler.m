function [t,y] = Euler(f,t,y0)

y = zeros(length(y0),length(t));
y(:,1) = y0;

for k=2:length(t)
    y(:,k) = y(:,k-1) + (t(k)-t(k-1))*f(t(k-1),y(:,k-1));
end

y = y';

end
    

