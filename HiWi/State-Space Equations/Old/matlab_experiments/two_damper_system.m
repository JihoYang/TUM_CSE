
clear
clc
close all
m1 = 10; % mass
k1 = 1; % spring constant
c1 = 1; % damping constant
m2 = 10;
k2 = 1;
c2 = 1;
Xin = @(t)(10*sin(10*t)); % External exiting displacement;
%F = @(t) (0);
Vin = @(t)(-100*cos(10*t));

x1(1) = 0;
v1(1) = 0;

x2(1) = 0;
v2(1) = 0;

dt = 0.0002;
t(1) = 0;

for n = 2:10000
    
    % For first mass 
    t(n) = (n-1)*dt;
    
    % Calculating Displacement
    x1(n) = x1(n-1) + dt*v1(n-1);
    
    % Calculating acceleration 
    acc1 = ( v1(n-1)*(c1+c2) + x1(n-1)*(k1+k2) - c2*v2(n-1) - k2*x2(n-1) - c1*Vin(t(n-1)) - k1*Xin(t(n-1)) )/m1;
    
    % Calculating Velocity
    v1(n) = v1(n-1) + dt*( acc1 );
    
    acc1(n) = acc1;
    f1(n) = m1*acc1(n); % Calculating Forces on first mass
    
    %
    % For second mass
    
    % Calculating Displacement
    x2(n) = x2(n-1) + dt*v2(n-1);
    
    % Calculating Acceleration
    acc2 = ( -c2*v2(n-1) -k2*x2(n-1) + k2*x1(n-1) + c2*v1(n-1) ) /m2;
    
    % Calculating Velocity
    v2(n) = v2(n-1) + dt*( acc2 );
    
    acc2(n) = acc2;
    f2(n) = m2*acc2(n); % Calculating Force on second mass
    
    
end


plot(t,x1);
hold on 
plot(t,x2,'k')

%figure(2);
%plot(t,acc2, 'r');
%figure(3);
%plot(t,v2, 'g');
%figure(4)
%plot(t,f2,'y');
