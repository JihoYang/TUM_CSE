% (5.4) Fehlerkriterium

ProtonTest
pause

% (5.5) Nicht-Negativität

RobertsonTest
pause

% (5.6) Steife AWP

ode15s(@proton,[0,8e5],[0;1;0]);
pause

ode45(@proton,[0,8e5],[0;1;0]);
pause

% steifes Testbeispiel

f = @(t,y) -1e15*y;
ode45(f,[0,1],1);
pause

ode23s(f,[0,1],1);
pause

% Beispiel mit Flammenzündung

delta = 1e-5; f = @(t,r) r^2-r^3;

ode45(f,[0,2/delta],delta);
pause

ode23s(f,[0,2/delta],delta);
pause

% 

% (5.7) Blow-Up

f = @(t,y) y^2;
ode45(f,[0,2],1);
pause

options = odeset('RelTol',1e-10);
ode45(f,[0,2],1,options);
pause