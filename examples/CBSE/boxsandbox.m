options.floating = true;
options.terrain = RigidBodyFlatTerrain();
options.dt = 0.01;
w = warning('off','Drake:TaylorVar:DoubleConversion');
r = TimeSteppingRigidBodyManipulator(fullfile(getDrakePath,'systems','plants','test','FallingBrickContactPoints.urdf'),options.dt,options);
feedback = FullStateFeedbackSensor();
r = addSensor(r, feedback);
body = 2;
xyz = [0; 0; 0];
rpy = zeros(3, 1); %nonzero rpy not supported in imu
imu = RigidBodyInertialMeasurementUnit(r, body, xyz, rpy);
r = r.addSensor(imu);
r = compile(r);

slamThing = SLAMThing(r, options);
sys = mimoCascade(r, slamThing);
options.floating = true;
options.terrain = RigidBodyFlatTerrain();
options.dt = 0.01;
w = warning('off','Drake:TaylorVar:DoubleConversion');
r = TimeSteppingRigidBodyManipulator(fullfile(getDrakePath,'systems','plants','test','FallingBrickContactPoints.urdf'),options.dt,options);
feedback = FullStateFeedbackSensor();
r = addSensor(r, feedback);
body = 2;
xyz = [0; 0; 0];
rpy = zeros(3, 1); %nonzero rpy not supported in imu
imu = RigidBodyInertialMeasurementUnit(r, body, xyz, rpy);
r = r.addSensor(imu);
r = compile(r);

slamThing = SLAMThing(r, options);
sys = mimoCascade(r, slamThing);



r = setSimulinkParam(r,'MinStep','0.001');
x0 = Point(r.getStateFrame);
x0.base_z = 10;
x0.base_zdot = 0.0;
x0.base_pitch = 0.0;
x0.base_roll = 0.5;
x0.base_yaw = 0.0;
x0.base_rolldot = 0.0;
x0.base_yawdot = 0.0;

v = r.constructVisualizer;
v.display_dt = .01;

if (1)
  % Run animation while it is simulating (as fast as possible, but probably
  % slower than realtime)
  s = warning('off','Drake:DrakeSystem:UnsupportedSampleTime');  % we are knowingly breaking out to a simulink model with the cascade on the following line.
  outs(1).system = 1;
  outs(1).output = 1;
  outs(2).system = 1;
  outs(2).output = 2;
  outs(3).system = 1;
  outs(3).output = 3;
  
  connection(1).from_output = 1;
  connection(1).to_input = 1;

  sys = mimoCascade(sys,v, connection, [], outs);
  warning(s);
  out = simulate(sys,[0 1],x0);
  
  xtraj_obj = DTTrajectory(out.tt, out.xx(1:12, :));
  xtraj_track = DTTrajectory(out.tt, out.xx(end-11:end, :));
  xtraj_obj = xtraj_obj.setOutputFrame(r.getStateFrame);
  xtraj_track = xtraj_track.setOutputFrame(r.getStateFrame);
  
  v.playback(xtraj_obj);
  v.playback(xtraj_track);
  
else 
% Run simulation, then play it back at realtime speed
  xtraj = simulate(r,[0 5],x0);
  v.playback(xtraj);
end


r = setSimulinkParam(r,'MinStep','0.001');
x0 = Point(r.getStateFrame);
x0.base_z = 10;
x0.base_zdot = 0.0;
x0.base_pitch = 0.0;
x0.base_roll = 0.5;
x0.base_yaw = 0.0;
x0.base_rolldot = 0.0;
x0.base_yawdot = 0.0;

v = r.constructVisualizer;
v.display_dt = .01;

if (1)
  % Run animation while it is simulating (as fast as possible, but probably
  % slower than realtime)
  s = warning('off','Drake:DrakeSystem:UnsupportedSampleTime');  % we are knowingly breaking out to a simulink model with the cascade on the following line.
  outs(1).system = 1;
  outs(1).output = 1;
  outs(2).system = 1;
  outs(2).output = 2;
  outs(3).system = 1;
  outs(3).output = 3;
  
  connection(1).from_output = 1;
  connection(1).to_input = 1;

  sys = mimoCascade(sys,v, connection, [], outs);
  warning(s);
  out = simulate(sys,[0 1],x0);
  
  xtraj_obj = DTTrajectory(out.tt, out.xx(1:12, :));
  xtraj_track = DTTrajectory(out.tt, out.xx(end-11:end, :));
  xtraj_obj = xtraj_obj.setOutputFrame(r.getStateFrame);
  xtraj_track = xtraj_track.setOutputFrame(r.getStateFrame);
  
  v.playback(xtraj_obj);
  v.playback(xtraj_track);
  
else 
% Run simulation, then play it back at realtime speed
  xtraj = simulate(r,[0 5],x0);
  v.playback(xtraj);
end