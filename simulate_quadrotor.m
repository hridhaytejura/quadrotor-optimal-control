
close all;
clear all;

time = 0;
%x0 = [2000, pi/4, -291, pi/16];
x0 = [100, pi/4, -10, pi/16];
m = 100;
I = 100;
g = 9.81;
L = 1;
c = 0.05;
Tend = 5;

A = [0 0 1 0;
     0 0 0 1;
     0 0 0 0;
     0 0 0 0];
B = [0 0;
     0 0;
     1/m 0;
     0 1/I];

q11 = 1;
q22 = 1;
q33 = 1;
q44 = 1;
Q = 100000* [q11 0 0 0;
     0 q22 0 0;
     0 0 q33 0;
     0 0 0 q44];
Q
% 100000*eye(4); % (assume Q is diagonal) x'Qx = q11 * x1^2 + q22 * x2^2 
% + q33 * x3^2 (velocity squared -> heat load squared)+ q44*x4^2
% more weight on state variable 

r11 = 1;
r22 = 1;

% Task: what is u'Ru equal to? 
R = eye(2); % Task: come up with a realistic cost model for fuel use

% u = -Kx (sometimes u = Kx)
K = lqr(A,B,Q,R);



% Ordinary Differential Equation (ODE) 4-5: 4th-5th order 
% 1st order: x(t+1) = x(t) + dx/dt * Delta_t ("1st derivative")
[t, x] = ode45(@(t,x) simple_quadrotor_dynamics(t, x, K, I, m, g, L, c), [0, Tend], x0)


figure();
hold on;
plot(t, x(:,1));
plot(t, (180/pi)*x(:,2));
plot(t, x(:,3));
plot(t, (180/pi)*x(:,4));
legend('y', '\theta', 'd/dt y', 'd/dt \theta');
title('Plot of Quadrotor states over time')