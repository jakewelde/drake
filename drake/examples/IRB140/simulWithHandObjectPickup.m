%% load in
fprintf('Setup...\n');
dt = 0.001;
options.dt = dt;
options.floating = false;
options.base_offset = [0;0;0]; %[-0.5, 0, 1.5]'; %For now positions on ground
options.base_rpy = [0, 0, 0]';
options.ignore_self_collisions = true;
options.collision = false;
options.hands = 'none';
r = IRB140('urdf/irb_140.urdf', options);
options_hand = options;
options_hand.hands = 'block_hand';
options_hand.terrain = RigidBodyFlatTerrain();
r_hand = IRB140('urdf/irb_140.urdf', options_hand);

pitch = 0.5;
yaw = 0.5;
rows = 100;
cols = 100;
rgbd_frame = RigidBodyFrame(1,[1;0;0.25],[0;0;pi],'rgbdframe');
r_hand = r_hand.addFrame(rgbd_frame);
rgbd_sensor = RigidBodyDepthSensor('rgbd', r_hand.findFrameId('rgbdframe'), -pitch, pitch, rows, -yaw, yaw, cols, 10.0);
rgbd_sensor = rgbd_sensor.enableLCMGL();
r_hand = r_hand.addSensor(FullStateFeedbackSensor());
r_hand = r_hand.addSensor(rgbd_sensor);
r_hand = r_hand.compile();

%hand_coll_links = {'palm'; 'finger_1_link_1'; 'finger_1_link_2'; 'finger_1_link_3';
%   'finger_2_link_1'; 'finger_2_link_2'; 'finger_2_link_3'; 'finger_middle_link_1'; 'finger_middle_link_2'; 'finger_middle_link_3'};
%box_link = 'drill_box';
%r_hand = r_hand.addLinksToCollisionFilterGroup(hand_coll_links, 'default', 2); 
%r_hand = r_hand.addLinksToCollisionFilterGroup(box_link, 'default', 3); 
%r_hand = r_hand.removeCollisionGroupsExcept({'default'});

options_cyl.floating = true;
r_hand = r_hand.addRobotFromURDF('drill_box.urdf', [0.4; 0.0; 0.05], [], options_cyl);

v=r.constructVisualizer();

%% initial config
% x0_hand = r_hand.getInitialState();
x0 = [-0.4;0.4;0.4;0;-0.8;0; zeros(6,1)];
x0_hand = r_hand.getInitialState();
x0_hand(1:12) = x0;

%% final pose
fprintf('Generating target traj\n');
options.visualize = true;

target_xtraj = runPlanning(x0(1:r.getNumPositions), [0.4, 0.5, 0.5]', options);


pd_control = irb140_trajfollow_block(r, target_xtraj);
clear ins; clear connection_1to2;
clear outs; clear connection_2to1;
%ins(1).system = 1;
%ins(1).input = 2;
for i=1:r_hand.getOutputFrame.getNumFrames
    outs(i).system = 1;
    outs(i).output = i;
end
connection_1to2(1).from_output = 1;
connection_1to2(1).to_input = 1;
connection_2to1(1).from_output = 1;
connection_2to1(1).to_input = 1;
sys = mimoFeedback(r_hand, pd_control, connection_1to2, connection_2to1, [],outs);

%% simulate
v=r_hand.constructVisualizer();
v.display_dt = 0.001;
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
output_select(2).system=1;
output_select(2).output=2;
sys = mimoCascade(sys,v,[],[],output_select);
warning(S);

traj = simulate(sys,[0 0.5],x0_hand);

% This doesn't see hand movements. Why?
playback(v,traj,struct('slider',true));