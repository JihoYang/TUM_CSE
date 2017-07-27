%---------------------------------------%
%-- This function contains test cases --%
%-- for newtons method                --%
%---------------------------------------%

close all;
clear al;
clc;

% Things that could make the method run slowly:
%   * Higher order roots
%   * Bad guess

% Things that could make he method diverge/not get a solution
%   * fd = 0
%   * Oscillating between two points 
%   * Bad guess

%TEST 1: Does it work at all?

f  = @(x) x^2 -1;
fd  = @(x) 2*x;
x0 = 3;
err = 1e-4;

%root = newtons_method(f,fd,x0,err)

%TEST 2: Stationary point: fd = 0 at a point.
x0 = 0 ;
%root = newtons_method(f,fd,x0,err)

x0 = pi/2;
f  = @(x) sin(x);
fd = @(x) cos(x);
%root = newtons_method(f,fd,x0,err)%zero is far off starting point

%TEST 4: Divergent
x0 = 2;
f  = @(x) x^(1/3);
fd = @(x) (1/3) * x^(-2/3); 
%root = newtons_method(f,fd,x0,err)

%TEST 5: No solution.
x0 = 1;
f  = @(x) x^2+1;
fd  = @(x) 2*x;
%root = newtons_method(f,fd,x0,err)

%TEST 6: No solution and oscillating
x0 = 1; 
f = @(x) exp(x) - 2*x;
fd = @(x) exp(x) - 2;
%root = newtons_method(f,fd,x0,err)

%TODO: %TEST 7: Poor initial guess.
%               A sin() which finds wrong root?

%TEST 8: High multiplicity
% Takes about 1300 steps for accuracy ~ 1e-4
x0 = 5;
f = @(x) (x-1)^10;
fd = @(x) 10*(x-1)^9;
root = newtons_method(f,fd,x0,err)

%TEST 9: Several roots

%TEST 10: Adaptive euler ???


