function sol=A29
vM=1; vC=2; eventTOL=1e-5; RelTol=1e-6; t0=0;
T=1; % "event" will happen before T

opt=odeset('RelTol',RelTol,'AbsTol',eventTOL*RelTol,...
    'Refine',1,'Events',@event);
sol=ode45(@rhs,[t0,T],[0;0;1;0],opt);

tt=unique([sol.x,linspace(t0,sol.xe,50)]); % points for plot

figure;
plot(sol.y(1,:),sol.y(2,:),'k*',deval(sol,tt,1),deval(sol,tt,2),'k',...
     sol.y(3,:),sol.y(4,:),'rx',deval(sol,tt,3),deval(sol,tt,4),'r');   
xlabel('x');ylabel('y');title('path of mouse and cat, resp.');
figure;
t=sol.x;
semilogy(t,sqrt(sum((deval(sol,t,[1,2])-deval(sol,t,[3,4])).^2)),'.');
xlabel('t');ylabel('distance');
 
function dz = rhs(~,z)
    % z(1:2) x/y position of mouse   z(3:4)   x/y position of cat
    r=z(1:2)-z(3:4);
    dz=[0;vM;vC*r/norm(r,2)];
end

function [value,isterminal,direction]=event(~,z)
    isterminal=1; direction=-1; 
    r=z(1:2)-z(3:4);
    value=norm(r,2)-eventTOL;
end

end