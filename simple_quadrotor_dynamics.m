% quadrotor: x = (y, theta, ydot, thetadot) 
% y is altitude
% theta is angle from vertical 
% ydot is vertical velocity
% thetadot is angular velocity
% x(1) = y, x(2) = theta, x(3) = ydot, x(4) = thetadot
% xdot = d/dt (x)

% F = (F1, F2)
% F1 = left propeller thrust
% F2 = right propeller thrust
% u(1) = total thrust = F1 + F2, u(2) = thrusting moment = L*(F2-F1)
% exercise: calculate F1, F2 from u1, u2

function xdot = simple_quadrotor_dynamics(t, x, K, I, m, g, L, c)
% t: time
% x: state
% K: control input gain 
% I: intertia (use 100)
% m: mass (use 100)
% g: gravity (use 9.81)
% L: thruster distance from center of mass (use 1)
% c: fluid resistance coefficient

u = -K*x + g

xdot(1,1) = x(3);
xdot(2,1) = x(4);
xdot(3,1) = u(1)/m * cos(x(2)) - g/m - c/m*x(3);
xdot(4,1) = u(2)/I;

end 
