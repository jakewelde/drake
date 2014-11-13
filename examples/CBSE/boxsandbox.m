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
  outs(1).output = 2;
  sys = mimoCascade(sys,v, [], [], outs);
  warning(s);
  out = simulate(sys,[0 2],x0);
else 
% Run simulation, then play it back at realtime speed
  xtraj = simulate(r,[0 5],x0);
  v.playback(xtraj);
end