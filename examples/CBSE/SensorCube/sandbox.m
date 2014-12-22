% Play with sensor cube data
load('lcmlog_2014_10_20_27.mat');

options.floating = true;
options.terrain = RigidBodyFlatTerrain();
options.dt = 0.01;
w = warning('off','Drake:TaylorVar:DoubleConversion');
r = TimeSteppingRigidBodyManipulator('SensorCubeContactPoints.urdf',options.dt,options);
feedback = FullStateFeedbackSensor();
r = addSensor(r, feedback);
body = 2;
xyz = [0; 0; 0];
rpy = zeros(3, 1); %nonzero rpy not supported in imu
imu = RigidBodyInertialMeasurementUnit(r, body, xyz, rpy);
r = r.addSensor(imu);
r = compile(r);

r = setSimulinkParam(r,'MinStep','0.001');
x0 = Point(r.getStateFrame);
x0.base_z = 1;
x0.base_zdot = 0.0;
x0.base_pitch = 0.0;
x0.base_roll = 0.5;
x0.base_yaw = 0.0;
x0.base_rolldot = 0.1;
x0.base_yawdot = 0.1;

v = r.constructVisualizer;
v.display_dt = .01;

% Run simulation, then play it back at realtime speed
%xtraj = simulate(r,[0 5],x0);
times = VICON_rlg_cbse_brick(:, 1).';
% fudge time points to make DTTrajectory happy -- it wants
% precisely regularly spaced time points, which we don't perfectly
% have, though we ain't fair from it
times = linspace(times(1), times(end), length(times));
times = (times - times(1))/1000000;
tstart = 7.5; tend = 9.0;
inds = times > tstart & times < tend;
times = times(inds);
poses = zeros(12, sum(inds));
poses(1:3, :) = VICON_rlg_cbse_brick(inds, 2:4).';
% set final pos at origin + half of block height
poses(1:3, :) = poses(1:3, :) - repmat(poses(1:3, end) - [0; 0; 0.05], [1, size(poses, 2)]);
% raise up a bit
poses(3, :) = poses(3, :) + 0.5;
indsi = find(inds);
imu_times = attitude(:, 1);
imu_times = (imu_times - imu_times(1))/1000000;
imu_inds = imu_times > tstart & imu_times < tend;
imu_data = attitude(imu_inds, :).';
for i=1:length(indsi)
  quat = VICON_rlg_cbse_brick(indsi(i), 5:8);
  poses(4:6, i) = quat2rpy(quat);
end
xtraj_constructed = DTTrajectory(times, poses);
xtraj_constructed = xtraj_constructed.setOutputFrame(v.getInputFrame);
v.playback(xtraj_constructed);

% so time range is
tf = tend - tstart;
% boundaries
x0 = xtraj_constructed.eval(xtraj_constructed.tt(1));
x1 = xtraj_constructed.eval(xtraj_constructed.tt(2));
xd0 = (x1 - x0) / (xtraj_constructed.tt(2) - xtraj_constructed.tt(1));
xfminus1 = xtraj_constructed.eval(xtraj_constructed.tt(end-1));
xf = xtraj_constructed.eval(xtraj_constructed.tt(end));
xdf = (xf - xfminus1) / (xtraj_constructed.tt(end) - xtraj_constructed.tt(end-1));
x0(7:end) = xd0(1:6);
xf(7:end) = xdf(1:6);

%% Simulate the block from the same init condition, but using our
% simpler, inelastic dynamics
xtraj_simul = simulate(r,[0 tend - tstart],x0);
v.playback(xtraj_simul);

% Now that we have it in a traj, go ahead and formulate as a
% trajectory optim problem
options = struct();
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

% # of samples we'll take
N=10;


scale_sequence = [10;1;0.5;0.1;.001;0];

for i=1:length(scale_sequence)
  scale = scale_sequence(i)
  
  options.compl_slack = scale*.01;
  options.lincompl_slack = scale*.001;
  options.jlcompl_slack = scale*.01;
  
  prog = ContactImplicitTrajectoryOptimization(r.getManipulator,N,tf,options);
  prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
  prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
  prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
  
  % final conditions constraint
  prog = addStateConstraint(prog,ConstantConstraint(x0(1:12)),1);
  %prog = addStateConstraint(prog,ConstantConstraint(xf(1:12)),N);
  
  % Add in IMU constraints
  break_inds = floor(linspace(1, sum(inds), N));
  break_inds_imu = floor(linspace(1, sum(imu_inds), N));
  for n=1:N
    % -> given J_imu at x_n, J_imu(x_n) * (x_n+1 - x_n) + J0 = z_n+1
    imu_constraint = zeros(6, 24);
    imu_equality = zeros(6, 1);
    x = poses(:, break_inds(n));
%     x_next = poses(:, break_inds(n+1));
     ximu_unsorted = imu_data(:, break_inds_imu(n));
%     ximu_unsorted_next = imu_data(:, break_inds_imu(n+1));
    ximu = zeros(10, 1);
    ximu(8:10) = ximu_unsorted(9:11);
    ximu(5:7) = ximu_unsorted(3:5);
%     ximu_next = zeros(10, 1);
%     ximu_next(8:10) = ximu_unsorted_next(9:11);
%     ximu_next(5:7) = ximu_unsorted_next(3:5);
    %     % Numerical:
    %     out_norm = imu.output(r.getManipulator,0,double(x), []);
    %     out_num = zeros(6, length(x));
    %     eps = 1E-6;
    %     add = x*0;
    %     for j=1:length(x)
    %       add = add*0;
    %       add(j) = eps;
    %       out_this = (imu.output(r.getManipulator,0,double(x+add), []) - out_norm)/eps;
    %       acc = out_this(8:10);
    %       gyro = out_this(5:7);
    %       out_num(:, j) = [acc; gyro];
    %     end
    %     out = out_num;
    %
    %     % 1st-order taylor expansion of imu output at x_n.
    %     % J * (xn+1 - x_n) + zn = zn+1
    %     % i.e. j*xn+1 = zn+1 - J0 + J*x_n
    %     % Jac * x_n v.playback(xtraj);+1 - Jac * x_n
    %     imu_constraint(:, 1:12) = -out;
    %     imu_constraint(:, 13:end) = out;
    %     % ... J0 = z_n+1
    %     imu_equality(1:3) = ximu_next(8:10) - ximu(8:10);
    %     imu_equality(4:6) = ximu_next(5:7) - ximu(5:7);
    %     delta = 0.1*scale;
    %     prog = prog.addStateConstraint(LinearConstraint(imu_equality-delta,imu_equality+delta,imu_constraint), {n:n+1});
    
    delta = 0.1*scale;
    prog = prog.addStateConstraint(FunctionHandleConstraint(ximu-delta,ximu+delta, 12, ...
                                  @(x_inp) imu.output(r.getManipulator,0,double(x_inp),[])), n);
    
  end
  
  clear traj_init;
  if i == 1
    traj_init.x = PPTrajectory(foh([0,tf],[x0(1:12),x0(1:12)]));
  else
    traj_init.x = xtraj;
    traj_init.l = ltraj;
  end
  tic
  [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
  toc
end
%%
v.playback(xtraj_constructed);
v.playback(xtraj);