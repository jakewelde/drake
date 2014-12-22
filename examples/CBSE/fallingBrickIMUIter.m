% uses mposa contact implicit trajectory optimization to generate traj
% for falling brick based on the final state.
rng(0)
visualize = true;

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
plant = RigidBodyManipulator(fullfile(getDrakePath,'systems','plants','test','FallingBrickContactPoints.urdf'),options);
warning(w);

x0 = [0;0;1.0; rand(3, 1); 0;0;4.9; rand(3, 1)];
N=10; tf=1.25;

plant_ts = TimeSteppingRigidBodyManipulator(plant,tf/(N-1));
feedback = FullStateFeedbackSensor();
plant_ts = addSensor(plant_ts, feedback);
body = 2;
xyz = [0; 0; 0];
rpy = zeros(3, 1); %nonzero rpy not supported in imu
%imu = RigidBodyInertialMeasurementUnit(plant_ts, body, xyz, rpy);
imu = TimeSteppingInertialMeasurementUnit(plant_ts, body, xyz, rpy);
plant_ts = plant_ts.addSensor(imu);
plant_ts = plant_ts.compile();
w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
x0WithSensor = [x0; x0; x0];
xtraj_ts = simulate(plant_ts,[0 tf],x0WithSensor);
x0 = xtraj_ts.eval(0);
xf = xtraj_ts.eval(xtraj_ts.tspan(end));
warning(w);
if visualize
  v = constructVisualizer(plant_ts);
  v.playback(xtraj_ts);
end

%%
imu_outs = [];
for i=1:N
  out = xtraj_ts.eval(xtraj_ts.tt(i));
  out = out(16:18);
  imu_outs = [imu_outs out];
end
figure;
plot(imu_outs.')
%%

options = struct();
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

scale_sequence = [0.1; .001; 0];
w = warning('off','Drake:FunctionHandleTrajectory');

total_eval_times = [];
for k=1:N-1
  fprintf('Using iter %d to %d\n', N-k, N);
  
  tic
  for i=1:length(scale_sequence)
    scale = scale_sequence(i)

    options.compl_slack = scale*.01;
    options.lincompl_slack = scale*.001;
    options.jlcompl_slack = scale*.01;

    ts = getBreaks(xtraj_ts);
    prog = ContactImplicitTrajectoryOptimization(plant,k+1,tf-ts(N-k),options);
    prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
    prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
    prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
    % prog = prog.setCheckGrad(true);

    %   snprint('snopt.out');

    % final conditions constraint
    %prog = addStateConstraint(prog,ConstantConstraint(x0(1:12)),1);
    prog = addStateConstraint(prog,ConstantConstraint(xf(1:12)),k+1);
    %prog = prog.addStateCost(FunctionHandleObjective(12, ...
    %        @(x_inp) norm(x_inp(1:12)-xf(1:12)), -1), N);

    % Add in IMU constraints
    dt = (xtraj_ts.tspan(end) - xtraj_ts.tspan(1))/(N-1);
    for n=N-k:N-1
      % -> given J_imu at x_n, J_imu(x_n) * (x_n+1 - x_n) + J0 = z_n+1
      imu_constraint = zeros(6, 24);
      imu_equality = zeros(6, 1);
      x = xtraj_ts.eval(ts(n));
      xlast = xtraj_ts.eval(ts(n+1));
      %xnext = xtraj_ts.eval(ts(n+1));
      ximu = x(13:end);
      %ximunext = xnext(13:enupdated);
      x = x(1:12);

      this_ind = n - (N-k-1);
      prog = prog.addStateCost(FunctionHandleObjective(24, ...
                                     @(x_inp_last, x_inp)  norm(imu.output(plant_ts,0,0,[double(x_inp); double(x_inp); double(x_inp_last)],[])-ximu), -1), ... 
                               [this_ind, this_ind + 1]);
    end

    clear traj_init;
    if (i == 1),
      if (k==1)
        traj_init.x = PPTrajectory(foh([0,tf-ts(N-k)],[xf(1:12),xf(1:12)]));
      else
        traj_init.x = PPTrajectory(foh([0, tf-ts(N-k+1)], [xtraj.eval(ts(N-k+1)), xtraj.eval(ts(N-k+1))]));
        traj_init.x = traj_init.x.append(xtraj.shiftTime(traj_init.x.tspan(end)));
      end
    else
      traj_init.x = xtraj;
      traj_init.l = ltraj;
    end
    %traj_init.x = xtraj_seed;
    w = warning('off','Drake:TaylorVar:DoubleConversion');
    
    [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf-ts(N-k),traj_init);
    
  end
  time = toc
  total_eval_times = [total_eval_times; k time];
end

if visualize
  %%
  v.playback(xtraj_ts);
  v.playback(xtraj);
end