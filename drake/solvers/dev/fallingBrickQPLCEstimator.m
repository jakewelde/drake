options.terrain = RigidBodyFlatTerrain();
options.floating = true;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
plant = RigidBodyManipulator(fullfile(getDrakePath,'systems','plants','test','FallingBrick.urdf'),options);
warning(w);
x0 = [0;0;.8;0.5; 0.1;pi/2;zeros(6,1)];

dt = 0.00333;
tf = 2;
plant_ts = TimeSteppingRigidBodyManipulator(plant,dt);
w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');

% Instrument it
pitch = 0.5;
yaw = 0.5;
rows = 100;
cols = 100;
scanOffset = [3;0;0.25;0;0;pi];
rgbd_frame = RigidBodyFrame(1,scanOffset(1:3),scanOffset(4:6),'rgbdframe');
plant_ts = plant_ts.addFrame(rgbd_frame);
rgbd_sensor = RigidBodyDepthSensor('rgbd', plant_ts.findFrameId('rgbdframe'), -pitch, pitch, rows, -yaw, yaw, cols, 10.0);
%rgbd_sensor = rgbd_sensor.enableLCMGL();
plant_ts = plant_ts.addSensor(FullStateFeedbackSensor());
plant_ts = plant_ts.addSensor(rgbd_sensor);
plant_ts = plant_ts.compile();

% Construct our QP estimator
%qplc_est = QPLCContactImplicitEstimator(plant_ts, x0, []);
%qplc_est = QPDART(plant_ts, x0, []);
qplc_est = QPMIContactImplicitEstimator(plant_ts, x0, []);
clear outs;
for i=1:plant_ts.getOutputFrame.getNumFrames
   outs(i).system = 1;
   outs(i).output = i;
end
for j=1:qplc_est.getOutputFrame.getNumFrames
   outs(i+j).system = 2;
   outs(i+j).output = j;
end
sys = mimoCascade(plant_ts, qplc_est, [], [], outs);

v = constructVisualizer(plant_ts);
% clear outs;
% for i=1:sys.getOutputFrame.getNumFrames
%     outs(i).system = 1;
%     outs(i).output = i;
% end
% sys = mimoCascade(sys, v, [], [], outs);
'starting sim'
ytraj = simulate(sys,[0 2],[x0; x0]);

gttraj = DTTrajectory(ytraj.tt,ytraj.xx(1:qplc_est.nx, :));
gttraj = gttraj.setOutputFrame(plant_ts.getOutputFrame.frame{1});
esttraj = DTTrajectory(ytraj.tt, ytraj.xx((qplc_est.nz+1):(qplc_est.nz+qplc_est.nx), :));
esttraj = esttraj.setOutputFrame(plant_ts.getOutputFrame.frame{1});

v.playback(gttraj);
v.playback(esttraj);