clear;
clc;
close all


% Stiffness values
k1 = 1;
k2 = 1;
k3 = 1;

% Damping coefficient values
c1 = 1;
c2 = 1;
c3 = 1;

% Masses of the bodies
m1 = 25.8;
m2 = 0;
m3 = 31.28;

% time Step
dt = 0.1;

% Initial positions of the masses
x1 = 7; y1 = 1; % M1 position
x2 = 5; y2 = 5; % M2 position
x3 = 2.5; y3 = 4.5; % M3 position

% Displacement evaluation Matrix
DEM = [1 0 0 dt 0 0 0 0 0;
       0 1 0 0 dt 0 0 0 0;
       0 0 1 0 0 dt 0 0 0];


% Velocity Evaluation Matrix
VEM = @(t)([-dt*(k1+k2)/m1 k1*dt/m1 k2*dt/m1 (1-(dt*(c1+c2)/m1)) dt*c1/m1 dt*c2/m1 dt*exp(-25*t)/m1 0 0;
     
        dt*k1/m2 -dt*(k1+k3)/m2 k3*dt/m2 dt*c1/m2 (1-(dt*(c1+c3)/m2)) dt*c3/m2 0 0 0;

        dt*k2/m3 dt*k3/m3 -dt*(k2+k3)/m3 dt*c2/m3 dt*c3/m3 (1-(dt*(c2+c3)/m3)) 0 0 0]);


% Setting up the initial displacement and velocity matrix
% First element in every column = x component
% Second component              = y component
% First column                  = mass 1 displacement
% Second column                 = mass 2 displacement
% third column                  = mass 3 displacement
% fourth column                 = mass 1 velocity
% fifth column                  = mass 2 velocity
% sixthh column                 = mass 3 velocity
% seventh column                = mass 1 force
% eight column                  = mass 2 force
% ninth column                  = mass 3 force
%
% Initially we move mass m1 with a displacement of 0.5 and see the
% vibrations.
dof = [];
dof(:,:,1) = [ 0 0 -1.5 0 0 0 0 0 0;
               0.5 0 0 0 0 0 0 0 0];
dof = dof';

% Time loop for 2000 time steps
figure;
for i = 2 : 5000
    
    % Calculating Displacements.
    displacement = DEM*dof(:,:,i-1);
    % Calculating Velocities.
    velocity = VEM(i*dt)*dof(:,:,i-1);
    % Calculating Forces.
    acceleration = velocity - dof(4:6,:,i-1);
    force = [dof(1,:,i-1)*m1/dt; dof(2,:,i-1)*m2/dt; dof(3,:,i-1)*m3/dt];
    % Assembling the degrees of freedom
    dof(:,:,i) = [displacement;velocity;force];
    dof(2,:,i) = [0;0];
    dof(5,:,i) = [0;0];
    %Recalculating the coordinates of the masses from displacemnts
    %calculated
    x1(i) = x1(1) + dof(1,1,i); y1(i) = y1(1) + dof(1,2,i);
    x2(i) = x2(1) + dof(2,1,i); y2(i) = y2(1) + dof(2,2,i);
    x3(i) = x3(1) + dof(3,1,i); y3(i) = y3(1) + dof(3,2,i);
    
   
        % Drawing circles to represent positions of the masses
        axis([1, 10, 0, 8]);
        title1 = sprintf('At Time :: %.3f',i*dt);
        title(title1);
        
        C1 = circle(x1(i-1),y1(i-1),0.20, 'r');  % Mass 1 in color red
        C2 = circle(x2(i-1),y2(i-1),0.1, 'b');  % Mass 2 in Blue color
        C3 = circle(x3(i-1),y3(i-1),0.1, 'k');  % Mass 3 in Black color
        
        % Plotting the force vectors
        x = [x1(i); x2(i); x3(i)];
        y = [y1(i); y2(i); y3(i)];
        x = [x;1];
        y = [y;1];
        fx = dof(4:6,1,i-1);
        fy = dof(4:6,2,i-1);
        %fx = acceleration(:,1);
        %fy = acceleration(:,2);
        fx = [fx;0.05];
        fy = [fy;0];
        Q = quiver(x,y,fx,fy);
        if(i == 2 || i == 3)
            pause();
        end
                
        pause(0.01);
        delete(C1);
        delete(C2);
        delete(C3);
        delete(Q);
    
    
end

