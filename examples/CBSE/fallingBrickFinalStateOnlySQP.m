% uses mposa contact implicit trajectory optimization to generate traj
% for falling brick based on the final state.
rng(0)
visualize = true;

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
plant = RigidBodyManipulator(fullfile(getDrakePath,'systems','plants','test','FallingBrickContactPoints.urdf'),options);
warning(w);
x0 = [0;0;1.0; rand(3, 1); 0;0;4.9; 2*rand(3, 1)];

N=10; tf=1.25;

plant_ts = TimeSteppingRigidBodyManipulator(plant,tf/(N-1));
plant_ts = plant_ts.compile();
w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
xtraj_ts = simulate(plant_ts,[0 tf],x0);
x0 = xtraj_ts.eval(xtraj_ts.tspan(1));
xf = xtraj_ts.eval(xtraj_ts.tspan(end));
warning(w);
if visualize
  v = constructVisualizer(plant_ts);
  v.playback(xtraj_ts);
end

options = struct();
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

scale_sequence = [1; 0.5; 0.1; 0.01; .001;0];
  tic
for i=1:length(scale_sequence)
  scale = scale_sequence(i)

  options.compl_slack = scale*.01;
  options.lincompl_slack = scale*.001;
  options.jlcompl_slack = scale*.01;
  
  prog = ContactImplicitTrajectoryOptimization(plant,N,tf,options);
  prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
  prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
  prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
  % prog = prog.setCheckGrad(true);
  
%   snprint('snopt.out');
  
  % final conditions constraint
 % prog = addStateConstraint(prog,ConstantConstraint(x0),1);
  prog = addStateConstraint(prog,ConstantConstraint(xf),N);
  %prog = prog.addStateCost(FunctionHandleObjective(12, ...
  %        @(x_inp) norm(x_inp-xf), -1), N);
  
  if i == 1,
    % start with weak guess
    traj_init.x = PPTrajectory(foh([0,tf],[xf,xf]));
  else
    traj_init.x = xtraj;
    traj_init.l = ltraj;
  end

  [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);

end
  toc
if visualize
  v.playback(xtraj);
end

% check if the two simulations did the same thing:
ts = getBreaks(xtraj_ts);
valuecheck(ts,getBreaks(xtraj));
xtraj_data = xtraj.eval(ts)
xtraj_ts_data = xtraj_ts.eval(ts)

