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
dt = 0.03; tf=2;

plant_ts = TimeSteppingRigidBodyManipulator(plant,dt);
feedback = FullStateFeedbackSensor();
plant_ts = addSensor(plant_ts, feedback);
body = 2;
xyz = [0; 0; 0];
rpy = zeros(3, 1); %nonzero rpy not supported in imu
%imu = RigidBodyInertialMeasurementUnit(plant_ts, body, xyz, rpy);
imu = TimeSteppingInertialMeasurementUnit(plant_ts, body, xyz, rpy, 0);
plant_ts = plant_ts.addSensor(imu);
plant_ts = plant_ts.compile();

% Add state estimator
options_se = options;
options_se.dt_state_est = dt;
options_se.state_frame_subind = 1;
full_state_est_handle = SharedDataHandle([]);
full_state_est_tt_handle = SharedDataHandle([]);
stateEstBlock = StateEstQPBlock(plant_ts, x0, options_se, full_state_est_tt_handle, full_state_est_handle);
sys = mimoCascade(plant_ts, stateEstBlock);

% Add a visualizer so we can watch
v = plant_ts.constructVisualizer;
v.display_dt = dt;
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

% kick it off
%w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
x0WithSensor = [x0; x0; x0];
out = simulate(sys,[0 tf],x0WithSensor);

xtraj_obj = DTTrajectory(out.tt, out.xx(1:12, :));
xtraj_track = DTTrajectory(out.tt, out.xx(end-11:end, :));
xtraj_track_final = DTTrajectory(full_state_est_tt_handle.getData(), full_state_est_handle.getData());
xtraj_obj = xtraj_obj.setOutputFrame(plant_ts.getStateFrame.getFrameByNum(1));
xtraj_track = xtraj_track.setOutputFrame(plant_ts.getStateFrame.getFrameByNum(1));
xtraj_track_final = xtraj_track_final.setOutputFrame(plant_ts.getStateFrame.getFrameByNum(1));

v.playback(xtraj_obj);
v.playback(xtraj_track);
v.playback(xtraj_track_final);