
    function [t,y]=G_plot(f,fd,y0,t0,dt,tmax)
    
    %Linearisation 1
    close all
    clc
    
    t = t0:dt:tmax; 
    dt = 0.25;
    y_new1 = [-20:20];
    y_new2 = y_new1;
    y_1 = y0*ones(size(y_new1));
    y_2 = y_1;
    
    
    g_1 = @ (y_new1) (1+7/2*dt.*(2-y_1/10)).*y_1./(1+7/20*dt*y_1) - y_new1;
     
     
     %Linearisation 2


    g_2 = @ (y_new2) (y_2+dt/2*(7*y_2.*(1-y_2/10)))./(1-7/2*dt*(1-y_2/10)) -y_new2;
    
    %Adams Moulton
     
    y = y0;
    g = @(x) x - y - dt/2 * (f(t,y) + f(t,x));
    
    
    figure (10);
    hold on;
    plot([-20 20], [0 0]);
    plot(y_new1, g_1(y_new1), 'r');
    plot(y_new2, g_2(y_new2), 'g');
    plot(y_new1, g(y_new1), 'b');

    
    end