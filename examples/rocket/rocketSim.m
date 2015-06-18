p = PlanarRocketPlant;
v = PlanarRocketVisualizer(p);
v.playback_speed = 2.0;
T = 10;
%u = setOutputFrame(PPTrajectory(zoh([0, T], [10 10; 0.1 0.1])), p.getInputFrame);
%sys = mimoCascade(u, p);
%xdes = [0;5;0;0;0;0];
%c = hoverLQR(p, xdes);

xtrajdes = setOutputFrame(PPTrajectory(foh([0, 3*T/4, T], [0 -10 -10; 15 5 5; 0 0 0; 0 0 0; 0 0 0; 0 0 0])), p.getOutputFrame);
c = hoverTVLQR(p, xtrajdes, T);
sys = feedback(p,c);

% trajopt = DircolTrajectoryOptimization(p,10,T);
% q_nom = p.getInitialState()*0;
% Q = ones(length(q_nom));
% b = zeros([length(q_nom), 1]);
% final_cost_function = QuadraticConstraint(-inf,inf,Q,b);
% trajopt = trajopt.addFinalCost(final_cost_function);
 
% [xtraj, utraj] = trajopt.solveTraj(T);

%sys = cascade(sys, v);

while (1)
  init = rand([6, 1])*2-1;
  init(1) = 0;
  init(2) = 15;
  init(3) = init(3)*0.1;
  init(4) = init(4)*2;
  init(5) = init(5)*2;
  init(6) = init(6)*0.1;
  
  sys.sys1 = sys.sys1.setInitialState(init);
  sys.sys1.noise = 1;
  xtraj = simulate(sys,[0 T]);
  v.playback(xtraj, struct('gui_control_interface', true));
  
  final_pos = xtraj.eval(T)
end