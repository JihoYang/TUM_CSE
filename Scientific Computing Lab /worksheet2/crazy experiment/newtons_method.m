%---------------------------------------%
%---------------------------------------%
%-- Root finding using Newtons method --%
%---------------------------------------%
%---------------------------------------%


% This implementation assumes that you have a continous derivative.
% NaN is returned if a point is not found.

% f: The function we want to find the roots of
% fd: The derivative of f
% x0: Initial guess for the root
% err: Tolerated error in the roots. err = |f(x) - f(approx_root) |
% root: The root found

function [root,lin1, lin2] = newtons_method(f,fd,x0,err, dt,lin1a, lin2a)

    % root_(1): Previous approx
    % root_(2): Current approx 
    % root_(3): Next approx
    % Initiating with different values to not trigger repeating check
    root_ = [x0+ 3*err x0 x0];


    step_limit = 1000;
    step = 0;

    y = [-20:20];
    y_1 = lin1a;
    y_2 = lin2a;
    g_1 = @ (y_new1) (1+7/2*dt.*(2-y_1/10)).*y_1./(1+7/20*dt*y_1) - y_new1;
    g_2 = @ (y_new2) (y_2+dt/2*(7*y_2.*(1-y_2/10)))./(1-7/2*dt*(1-y_2/10)) -y_new2;

    lin1 = (1+7/2*dt.*(2-y_1/10)).*y_1./(1+7/20*dt*y_1);
    lin2 = (y_2+dt/2*(7*y_2.*(1-y_2/10)))./(1-7/2*dt*(1-y_2/10));
    
    figure;
    hold on;
    while(abs(f(root_(3))) > err && f(root_(3)) ~= 0)
        step = step + 1;

        % Stops if the method uses too many steps
        if(step > step_limit)
            %disp('Too many steps');
            %disp('No solution or divergent?');
            root = NaN;
            return;
        end

        
        % Test if fd = 0,Inf or -Inf
        %        * zero is defined as below our accuracy limit.
        if(abs(fd(root_(3))) < err ...
           || fd(root_(3)) == Inf || fd(root_(3)) == -Inf)
            %disp('fd = 0');

            %Should add something on the relative scale of root_(3) 
            %instead of just err.
            root_(3) = root_(3) + err;
        end


        % Test for repeating values. 
        if( root_(1) == root_(3))
            %disp('Repeating values');
            root_(3) = (root_(3) + root_(1))/2;

        end %if

        %Calculates the root_ values for next step;
        root_(1) = root_(2);
        root_(2) = root_(3);
        root_(3) = root_(3)-f(root_(3))/fd(root_(3));
        

        
        
        pause(0.05);
        
        
        
    end %while
    plot(y, f(y), 'r');
    plot(y, g_1(y), 'g');
    plot(y, g_2(y), 'b');
        
    box on;
    grid on;
    root = root_(3);
    step;
    pause;
    hold off;
end %function


