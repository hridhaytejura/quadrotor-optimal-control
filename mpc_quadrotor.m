clear all
close all

time = 0;
x0 = [100, pi/4, -10, pi/16];
x0 = [2000, pi/4, -291, pi/16];
m = 100;
I = 100;
g = 9.81;
L = 1;
c = 0.05;
Tend = 10;

A = [0 0 1 0;
     0 0 0 1;
     0 0 0 0;
     0 0 0 0];
B = [0 0;
     0 0;
     1/m 0;
     0 1/I];

% LQR and MPC weights
q11 = 1;
q22 = 1;
q33 = 1;
q44 = 1;
Q = 100000* [q11 0 0 0;
     0 q22 0 0;
     0 0 q33 0;
     0 0 0 q44];



r11 = 1;
r22 = 1;

% Task: what is u'Ru equal to? 
R = [r11 0 ;
    0 r22]; % Task: come up with a realistic cost model for fuel use

%%
% MPC Formulation

% timing parameters
Ts = .01;
Tend = 5;
tspan = 0:Ts:Tend;
N = length(tspan);
Tlook = 1;
Nlook = floor(Tlook/Ts);

% discretize
[Ad,Bd] = c2d(A,B,Ts);

% set initial conditions and target r


x0 = [100, pi/4, -10, pi/16]';
xref = [0, 0, 0, 0]';
r = xref;
dx0 = x0 - xref;

rhat = repmat(r,Nlook,1);
xrefhat = repmat(xref, Nlook, 1);

% initialize
uhat = zeros(Nlook,1);
u = 0;
xk = x0;
dxk = dx0;

% use P=dare(...) for constrained LQR formulation
% P = Q;
P = idare(Ad,Bd,Q,R);
% P = icare(A,B,Q,R) / Ts;
%
% Iterative method for solving DARE
% P = ones(4,4);
% Plist = [];
% for i = 1:1000
%   P = Q + Ad'*P*Ad - Ad'*P*Bd*inv(R+Bd'*P*Bd)*Bd'*P*Ad;
%   if mod(i,100) == 1
%     disp(P)
%   end    
%   Plist(:,i) = P(:);
% end
% plot(Plist')

% expand cost matrices
Qhat = {};
Rhat = {};
for i = 1:Nlook
  Qhat{i} = Q;
  Rhat{i} = R;
  if i == Nlook
    Qhat{i} = P;
  end
end
Qhat = blkdiag(Qhat{:});
Rhat = blkdiag(Rhat{:});

% Prediction model
nx = size(Ad,1);
nu = size(Bd,2);
E = [];
F  = [];
Apow  = eye(nx);
F_row = zeros(nx, Nlook*nu);
for i = 1:Nlook
  idx = (nu*(i-1)+1):(nu*i);
  F_row = Ad * F_row;
  F_row(:,idx) = Bd;
  F = [F; F_row];

  Apow = Ad*Apow;
  E = [E; Apow];
end

% State constraints: wall 
ground = -10; % -inf, -10, -1
maxrotrate = 10;
xmin = [ground -inf -inf -inf]';  
xmax = [ inf  inf inf  maxrotrate]';
xmax = repmat(xmax, Nlook, 1);
xmin = repmat(xmin, Nlook, 1);
Ix = [eye(Nlook*nx); -eye(Nlook*nx)];
Aineq_x = Ix*F;
bineq_x = [xmax; -xmin] - Ix * (E*dxk + xrefhat); 

% Input constraints
umin = -inf; % -inf, -40, -100;
umax = inf;  % inf, 40, 10;
Iu = [eye(Nlook*nu); -eye(Nlook*nu)];
Aineq_u = Iu;
bineq_u = [umax*ones(Nlook,1); -umin*ones(Nlook,1)];
  

%%
quadprog_opts =  optimset('Display','off');
iplot = 0;
cmap = lines;

% Simulate MPC
for i = 1:N
  
  fprintf('Step %d of %d \n', i, N)
 
  H = F'*Qhat*F + Rhat;
  H = (H+H')/2;
  ft = (dxk'*E' + (xrefhat-rhat)') * Qhat * F;
  f = ft';  
  
  % inequality constraints
  bineq_x = [xmax; -xmin] - Ix * (E*dxk + xrefhat);
  Aineq = [Aineq_x; Aineq_u];
  bineq = [bineq_x; bineq_u];
  iuse = find(~isinf(bineq));
  bineq = bineq(iuse);
  Aineq = Aineq(iuse,:);
  
  
  % Solve for solution
  % uhat = -inv(H)*f;  % Unconstrained solution
  uhat = quadprog(H,f,Aineq,bineq,[],[],[],[],[],quadprog_opts);
  
  
  u = uhat(1);
  
  % simulate
  t_interval = 0:.001:Ts;
  [t,x] = ode45(@(t,x) simple_quadrotor_dynamics(t, x, u, I, m, g, L, c), t_interval, xk)
  xk = x(end,:)';
  dxk = xk - xref;
  
  xmpc(i,:) = xk;
  umpc(i) = u;
  
   
  % Plot predictions along the way
  if mod(i,20) == 1
    iplot = iplot + 1;
    
    umpc_predicted = uhat;
    umpc_predicted = reshape(umpc_predicted, nu, []);

    xmpc_predicted = E*dxk + F*uhat + xrefhat;
    xmpc_predicted = reshape(xmpc_predicted, nx, []);
    tmpc_predicted = tspan(i) + (Ts:Ts:Tlook);
    
    figure(1)
    c = cmap(iplot, :);
    subplot(211)
    hold on
    scatter(tspan(i), u, 'markerfacecolor', c);
    plot(tmpc_predicted, umpc_predicted, '--', 'color', c) 
    plot(tspan(1:i), umpc, 'k')
    
    subplot(212)
    hold on
    scatter(tspan(i) * ones(4,1), xk, 'markerfacecolor', c);
    plot(tmpc_predicted, xmpc_predicted, '--', 'color', c) 
    plot(tspan(1:i), xmpc, 'k')
    
    if i == 1
      umpc_predicted1 = umpc_predicted;
      xmpc_predicted1 = xmpc_predicted;
      tmpc_predicted1 = tmpc_predicted;
    
      figure(1)
      subplot(211)
      grid on
      plot(tmpc_predicted, umpc_predicted, 'k')
      subplot(212)
      grid on
      plot(tmpc_predicted, xmpc_predicted, 'k')
    end
  end
end


%%
% Simulate LQR
K = lqr(A,B,Q,R);

uclip = @(u) min(max(u, umin),umax);
%[t,xlqr] = ode45(@(t,x)cartpend(x,m,M,L,g,d,uclip(-K*(x-r))),tspan,x0);
[t, xlqr] = ode45(@(t,x) simple_quadrotor_dynamics(t, x, uclip(-K*(x-r)), I, m, g, L, c), tspan, x0)
ulqr = uclip(-K*(xlqr'-r));

%% 
% Plot results

% figure
% subplot(211)
% grid on
% title('u', 'fontsize', 18)
% hold on
% plot(tspan, umpc, 'b', 'linewidth', 2)
% plot(tspan, ulqr, 'c', 'linewidth', 2)
% plot(tmpc_predicted1, umpc_predicted1, '--r', 'linewidth', 2)
% if ~isinf(umin)
%   yline(umin, 'r', 'linewidth', 2)
%   yline(umax, 'r', 'linewidth', 2)
% end
% legend('MPC', 'LQR', 'PREDICTION @ t=0', 'fontsize', 16)
% % ylim([-40 40])
% xlim([tspan(1) tspan(end)])
% 
% subplot(212)
% grid on
% title('x', 'fontsize', 18)
% hold on
% plot(tspan, xmpc, 'b', 'linewidth', 2)
% plot(tspan, xlqr, 'c', 'linewidth', 2)
% plot(tmpc_predicted1, xmpc_predicted1', '--r', 'linewidth', 2)
% xlim([tspan(1) tspan(end)])
% if ~isinf(wall)
%   yline(wall, 'r', 'linewidth', 2)
% end
% xlabel('Time [s]')


%%
% figure
% klist = floor(linspace(1, N, 100));

% for k = klist

%   drawcartpend_bw(xmpc(k,:),m,M,L,wall);
% %   drawcartpend_bw(xlqr(k,:),m,M,L,wall);
%   
%   %   filename = 'cartpend_mpc.gif';
%   %   frame = getframe(gcf);
%   %   im = frame2im(frame);
%   %   [imind,cm] = rgb2ind(im,256);
%   %   if k == 1
%   %       imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.02);
%   %   else
%   %       imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.02);
%   %   end
% end


% function drawcartpend_bw(y,m,M,L,wall_x)
%   x = y(1);
%   th = y(3);

  % kinematics
  % x = 3;        % cart position
  % th = 3*pi/2;   % pendulum angle

  % dimensions
  % L = 2;  % pendulum length
%   W = 1*sqrt(M/5);  % cart width
%   H = .5*sqrt(M/5); % cart height
%   wr = .2; % wheel radius
%   mr = .3*sqrt(m); % mass radius
% 
%   % positions
%   % y = wr/2; % cart vertical position
%   y = wr/2+H/2; % cart vertical position
%   w1x = x-.9*W/2;
%   w1y = 0;
%   w2x = x+.9*W/2-wr;
%   w2y = 0;
% 
%   px = x + L*sin(th);
%   py = y - L*cos(th);
% 
%   plot([-10 10],[0 0],'w','LineWidth',2)
%   hold on
%   rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1])
%   rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
%   rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
% 
%   % rectangle('Position',
% 
%   plot([x px],[y py],'w','LineWidth',2)
% 
%   rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.3 0.3 1],'EdgeColor',[1 1 1])
% 
%   % WALL CONSTRAINT
%   if exist('wall_x', 'var') && ~isinf(wall_x)
%     xline(wall_x - 0.05 - W/2, 'color', 'w', 'Linewidth', 4)
%   end
% 
%   % set(gca,'YTick',[])
%   % set(gca,'XTick',[])
%   xlim([-5 5]);
%   ylim([-2 2.5]);
%   set(gca,'Color','k','XColor','w','YColor','w')
%   set(gcf,'Position',[10 900 800 400])
%   set(gcf,'Color','k')
%   set(gcf,'InvertHardcopy','off')
% 
%   % box off
%   drawnow
%   hold off
% end




















