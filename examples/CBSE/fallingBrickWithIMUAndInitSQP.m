% uses mposa contact implicit trajectory optimization to generate traj
% for falling brick based on the final state.
rng(0)
visualize = true;

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
plant = RigidBodyManipulator(fullfile(getDrakePath,'systems','plants','test','FallingBrickContactPoints.urdf'),options);
warning(w);
x0 = [0;0;2.0;0.1*randn(3,1); 0; 0; -2.0; 1; 0; 0];

N=10; tf=2.0;

plant_ts = TimeSteppingRigidBodyManipulator(plant,tf/(N-1));
feedback = FullStateFeedbackSensor();
plant_ts = addSensor(plant_ts, feedback);
body = 2;
xyz = [0; 0; 0];
rpy = zeros(3, 1); %nonzero rpy not supported in imu
imu = RigidBodyInertialMeasurementUnit(plant_ts, body, xyz, rpy);
plant_ts = plant_ts.addSensor(imu);
plant_ts = plant_ts.compile();
w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
xtraj_ts = simulate(plant_ts,[0 tf],x0);
x0 = xtraj_ts.eval(0);
xf = xtraj_ts.eval(xtraj_ts.tspan(end));
warning(w);
if visualize
  v = constructVisualizer(plant_ts);
  v.playback(xtraj_ts);
end

options = struct();
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

scale_sequence = [1;.001;0];

for i=1:length(scale_sequence)
  scale = scale_sequence(i);

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
  prog = addStateConstraint(prog,ConstantConstraint(x0(1:12)),1);
  %prog = addStateConstraint(prog,ConstantConstraint(xf(1:12)),N);
%   
  % Add in IMU constraints
  ts = getBreaks(xtraj_ts);
  dt = (xtraj_ts.tspan(end) - xtraj_ts.tspan(1))/N;
  for n=1:N
    % -> given J_imu at x_n, J_imu(x_n) * (x_n+1 - x_n) + J0 = z_n+1
    imu_constraint = zeros(6, 24);
    imu_equality = zeros(6, 1);
    x = xtraj_ts.eval(ts(n));
    ximu = x(13:end);
    x = x(1:12);
    % Numerical:
    out_norm = imu.output(plant_ts.getManipulator,0,double(x), []);
    out_num = zeros(6, length(x));
    eps = 1E-6;
    add = x*0;
    for j=1:length(x)
      add = add*0;
      add(j) = eps;
      out_this = (imu.output(plant_ts.getManipulator,0,double(x+add), []) - out_norm)/eps;
      acc = out_this(8:10);
      gyro   % Add in IMU constraints
  ts = getBreaks(xtraj_ts);
  dt = (xtraj_ts.tspan(end) - xtraj_ts.tspan(1))/N;
  for n=1:N
    % -> given J_imu at x_n, J_imu(x_n) * (x_n+1 - x_n) + J0 = z_n+1
    imu_constraint = zeros(6, 24);
    imu_equality = zeros(6, 1);
    x = xtraj_ts.eval(ts(n));
    ximu = x(13:end);
    x = x(1:12);
    % Numerical:
    out_norm = imu.output(plant_ts.getManipulator,0,double(x), []);
    out_num = zeros(6, length(x));
    eps = 1E-6;
    add = x*0;
    for j=1:length(x)
      add = add*0;
      add(j) = eps;
      out_this = (imu.output(plant_ts.getManipulator,0,double(x+add), []) - out_norm)/eps;
      acc = out_this(8:10);
      gyro = out_this(5:7);
      out_num(:, j) = [acc; gyro];
    end
    acc = out_norm(8:10);
    gyro = out_norm(5:7);
    out = out_num;

    % 1st-order taylor expansion of imu output at x_n.
    % J * (xn+1 - x_n) + zn = zn+1
    % i.e. j*xn+1 = zn+1 - J0 + J*x_n
    % Jac * x_n+1 - Jac * x_n
    %if (n~=1)
      imu_constraint(:, 1:12) = -out;
      imu_constraint(:, 13:end) = out;
   % end
    %imu_constraint(offset_inner+1:offset_inner+6, offset+1:offset+12) = -out;
    % ... J0 = z_n+1
    %imu_constraint(offset_inner+1:offset_inner+6, offset+13:offset+18) = eye(6);
    imu_equality(1:3) = (ximu(8:10) - acc);
    imu_equality(4:6) = (ximu(5:7) - gyro);
    xind = (n-1)*12+1:(n+1)*12;
    prog = prog.addConstraint(LinearConstraint(imu_equality,imu_equality,imu_constraint), xind);
 end
= out_this(5:7);
      out_num(:, j) = [acc; gyro];
    end
    acc = out_norm(8:10);
    gyro = out_norm(5:7);
    out = out_num;

    % 1st-order taylor expansion of imu output at x_n.
    % J * (xn+1 - x_n) + zn = zn+1
    % i.e. j*xn+1 = zn+1 - J0 + J*x_n
    % Jac * x_n+1 - Jac * x_n
    %if (n~=1)
      imu_constraint(:, 1:12) = -out;
      imu_constraint(:, 13:end) = out;
   % end
    %imu_constraint(offset_inner+1:offset_inner+6, offset+1:offset+12) = -out;
    % ... J0 = z_n+1
    %imu_constraint(offset_inner+1:offset_inner+6, offset+13:offset+18) = eye(6);
    imu_equality(1:3) = (ximu(8:10) - acc);
    imu_equality(4:6) = (ximu(5:7) - gyro);
    xind = (n-1)*12+1:(n+1)*12;
    prog = prog.addConstraint(LinearConstraint(imu_equality,imu_equality,imu_constraint), xind);
 end
  
  
  if i == 1,
    traj_init.x = PPTrajectory(foh([0,tf],[xf(1:12),xf(1:12)]));
  else
    traj_init.x = xtraj;
    traj_init.l = ltraj;
  end
  tic
  [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
  toc
end

if visualize
  v.playback(xtraj);
end

% check if the two simulations did the same thing:
ts = getBreaks(xtraj_ts);
valuecheck(ts,getBreaks(xtraj));
xtraj_data = xtraj.eval(ts)
xtraj_ts_data = xtraj_ts.eval(ts)

