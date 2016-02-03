options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.multiple_contacts = false;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
plant = PlanarRigidBodyManipulator(fullfile(getDrakePath,'systems','plants','test','FallingBrickContactPoints.urdf'),options);
warning(w);
x0 = [0;1;0; 0; 0; 0.2];

dt = 0.00333;
tf = 2;
plant_ts = TimeSteppingRigidBodyManipulator(plant,dt);
w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');

% Instrument it
pitch = 0.5;
yaw = 0.0;
rows = 100;
cols = 1;
scanOffset = [3;0;0.25;0;0;pi];
rgbd_frame = RigidBodyFrame(1,scanOffset(1:3),scanOffset(4:6),'rgbdframe');
plant_ts = plant_ts.addFrame(rgbd_frame);
rgbd_sensor = RigidBodyDepthSensor('rgbd', plant_ts.findFrameId('rgbdframe'), -pitch, pitch, rows, -yaw, yaw, cols, 10.0);
plant_ts = plant_ts.addSensor(FullStateFeedbackSensor());
plant_ts = plant_ts.addSensor(rgbd_sensor);
plant_ts = plant_ts.compile();

PF = PFModeHopping(plant_ts, 10, x0);
clear outs;
for i=1:plant_ts.getOutputFrame.getNumFrames
   outs(i).system = 1;
   outs(i).output = i;
end
for j=1:PF.getOutputFrame.getNumFrames
   outs(i+j).system = 2;
   outs(i+j).output = j;
end
sys = mimoCascade(plant_ts, PF, [], [], outs);

v = constructVisualizer(plant_ts);
v.xlim = [-2 2]; v.ylim = [-0.5 3];
clear cons;
clear outs;
cons(1).from_output = 1;
cons(1).to_input = 1;
for i=1:sys.getOutputFrame.getNumFrames
 outs(i).system = 1;
 outs(i).output = i;
end
sys = mimoCascade(sys, v, cons, [], outs);
'starting sim'
ytraj = simulate(sys,[0 2],[x0; PF.getInitialState()]);