options.floating = true;
options.terrain = RigidBodyFlatTerrain();
% w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
% warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
% warning('off','Drake:RigidBodyManipulator:BodyHasZeroInertia');
% warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
% warning('off','Drake:RigidBodyGeometry:SimplifiedCollisionGeometry');
% warning('off','Drake:RigidBodyManipulator:ReplacedCylinder');
r = PlanarRigidBodyManipulator('box2d.urdf',options);
% warning(w);

%r = setSimulinkParam(r,'MinStep','0.001');
x0 = Point(r.getStateFrame);
x0.base_z = 1;
x0.base_zdot = -1;
x0.base_relative_pitch = 0.5;

v = r.constructVisualizer;
v.display_dt = .001;

if (1)
  % Run animation while it is simulating (as fast as possible, but probably
  % slower than realtime)
  s = warning('off','Drake:DrakeSystem:UnsupportedSampleTime');  % we are knowingly breaking out to a simulink model with the cascade on the following line.
  sys = cascade(r,v);
  warning(s);
  simulate(sys,[0 10],x0);
else 
% Run simulation, then play it back at realtime speed
  xtraj = simulate(r,[0 5],x0);
  v.playback(xtraj);
end