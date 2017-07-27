clear
clc
close all
m = 10; % massclear
k = 5; % spring constant
c = 5; % damping constant
%F = @(t)(-sin(10*t)); % External exiting force;
F = @(t) (0);
x(1) = -5;
v(1) = 0;
force = 0;
dt = 0.005;
t(1) = 0;

% 100 time steps

for n = 2:5001
    
    t(n) = (n-1)*dt;
    
    if(t(n)<0.2)
           f = 0;
    else
        f = F(n-1*dt);
    end
    
    x(n) = x(n-1) + dt*v(n-1);
    dummy = (f - c*v(n-1) - k*x(n-1))/m ;
    v(n) = v(n-1) + dt*( dummy );
    acc(n) = dummy;
    
    force(n) = m*acc(n);
  
end

plot(t,x);
hold on
plot(t,acc, 'r');
plot(t,v, 'g');
plot(t,force,'k');


