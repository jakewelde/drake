function rocketSim

p = PlanarRocketPlant;
v = PlanarRocketVisualizer(p);
ground = 5;
v = v.addPlatform(0, ground-p.ldraw/2, 2);
v.playback_speed = 2.0;
T = 10;

%u = setOutputFrame(PPTrajectory(zoh([0, T], [10 10; 0.1 0.1])), p.getInputFrame);
%sys = mimoCascade(u, p);
%xdes = [0;5;0;0;0;0];
%c = hoverLQR(p, xdes);

%xtrajdes = setOutputFrame(PPTrajectory(foh([0, 3*T/4, T], [0 -10 -10; 15 5 5; 0 0 0; 0 0 0; 0 0 0; 0 0 0])), p.getOutputFrame);
%c = hoverTVLQR(p, xtrajdes, T);
%sys = feedback(p,c);

N = 20;
trajopt = DircolTrajectoryOptimization(p,N,T);
%trajopt = trajopt.addFinalCost(@finalCost);
trajopt = trajopt.addRunningCost(@cost);
x0 = [-30;50;0;0;0;0];
x0(1) = x0(1) + 5*(rand*2-1);
x0(4:6) = x0(4:6) + rand(3, 1)*2-1;
xdes = [0;ground;0;0;0;0];
trajopt = trajopt.addStateConstraint(ConstantConstraint(double(x0)),1);
trajopt = trajopt.addStateConstraint(ConstantConstraint(double(xdes)),N);
for i=1:N
  % input limits
  trajopt = trajopt.addInputConstraint(BoundingBoxConstraint([0;-0.2], [30; 0.2]), i);
  % must always be above target
  trajopt = trajopt.addStateConstraint(BoundingBoxConstraint([-inf; ground; -inf; -inf; -inf; -inf], ...
    [inf; inf; inf; inf; inf; inf]), i);
  % Horizontal velocity must be low near the ground
  alpha = 0.05; % rads allowed per meter above "ground"
  A = [0 alpha -1 0 0 0;
       0 -alpha -1 0 0 0;];
  lb = [alpha*ground; -inf]; ub = [inf; -alpha*ground];
  trajopt = trajopt.addStateConstraint(LinearConstraint(lb, ub, A), i);
end
[xtraj, utraj] = trajopt.solveTraj(T);
xtraj = xtraj.setOutputFrame(v.getInputFrame);
v.playback(xtraj, struct('gui_control_interface', true));

% now stabilize with tvlqr
c = hoverTVLQR(p, xtraj, T);
sys = feedback(p,c);
sys.sys1 = sys.sys1.setInitialState(x0);
sys.sys1.noise = 1;
xtraj = simulate(sys,[0 T]);
v.playback(xtraj, struct('gui_control_interface', true));

%sys = cascade(sys, v);
% 
% while (1)
%   init = rand([6, 1])*2-1;
%   init(1) = 0;
%   init(2) = 15;
%   init(3) = init(3)*0.1;
%   init(4) = init(4)*2;
%   init(5) = init(5)*2;
%   init(6) = init(6)*0.1;
%   
%   sys.sys1 = sys.sys1.setInitialState(init);
%   sys.sys1.noise = 1;
%   xtraj = simulate(sys,[0 T]);
%   v.playback(xtraj, struct('gui_control_interface', true));
%   
%   final_pos = xtraj.eval(T)
% end

end

function [g,dg] = cost(dt,x,u)
Q = diag([0, 0, 10.0, 100.0, 100.0, 100.0]);
R = [1000, 0.0; 0.0, 1000];
g = x'*Q*x + u'*R*u;
%g = sum((R*u).*u,1);
dg = [0, 2*x'*Q,2*u'*R];
%dg = zeros(1, 1 + size(x,1) + size(u,1));
end