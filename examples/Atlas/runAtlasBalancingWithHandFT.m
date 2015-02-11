function runAtlasBalancingWithHandFT(use_mex)

if ~checkDependency('gurobi')
  warning('Must have gurobi installed to run this example');
  return;
end

path_handle = addpathTemporary(fullfile(getDrakePath(), 'examples', 'ZMP'));
import atlasControllers.*;

% put robot in a random x,y,yaw position and balance for 2 seconds
visualize = true;

if (nargin<1); use_mex = true; end
if (nargin<2); use_angular_momentum = false; end

% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

options.floating = true;
options.dt = 0.002;
options.num_external_force_frames = 1;
r = Atlas('urdf/atlas_minimal_contact.urdf',options);
options_hand = options;
options_hand.hands = 'robotiq_weight_only'; % so the actual system being controlled will have
% an extra kg or two on the right hand
r_actual = Atlas('urdf/atlas_minimal_contact.urdf', options_hand);
r = r.removeCollisionGroupsExcept({'heel','toe'});
r_actual = r_actual.removeCollisionGroupsExcept({'heel', 'toe'});
r = compile(r);
r_actual = compile(r_actual);

nq = getNumPositions(r);

% set initial state to fixed point
load('data/atlas_fp.mat');
xstar(1) = 0.1*randn();
xstar(2) = 0.1*randn();
xstar(6) = pi*randn();
r = r.setInitialState(xstar);

x0 = xstar;
q0 = x0(1:nq);
kinsol = doKinematics(r,q0);

com = getCOM(r,kinsol);

% build TI-ZMP controller
footidx = [findLinkId(r,'r_foot'), findLinkId(r,'l_foot')];
foot_pos = terrainContactPositions(r,kinsol,footidx);
comgoal = mean([mean(foot_pos(1:2,1:4)');mean(foot_pos(1:2,5:8)')])';
limp = LinearInvertedPendulum(com(3));
[~,V] = lqr(limp,comgoal);

foot_support = RigidBodySupportState(r,find(~cellfun(@isempty,strfind(r.getLinkNames(),'foot'))));


ctrl_data = QPControllerData(false,struct(...
  'acceleration_input_frame',atlasFrames.AtlasCoordinates(r),...
  'D',-com(3)/9.81*eye(2),...
  'Qy',eye(2),...
  'S',V.S,...
  's1',zeros(4,1),...
  's2',0,...
  'x0',[comgoal;0;0],...
  'u0',zeros(2,1),...
  'y0',comgoal,...
  'qtraj',x0(1:nq),...
  'support_times',0,...
  'supports',foot_support,...
  'mu',1.0,...
  'ignore_terrain',false,...
  'constrained_dofs',[]));

% instantiate QP controller
options.slack_limit = 30.0;
options.w_qdd = 0.001*ones(nq,1);
options.w_grf = 0;
options.w_slack = 0.001;
options.debug = false;
options.use_mex = use_mex;

if use_angular_momentum
  options.Kp_ang = 1.0; % angular momentum proportunal feedback gain
  options.W_kdot = 1e-5*eye(3); % angular momentum weight
else
  options.W_kdot = zeros(3);
end

options.use_mex = 0;
% add the hand ft frame to right hand
r_hand_body = findLinkId(r,'r_hand');
qp = QPController(r,{},{atlasFrames.ExternalForceTorque()},ctrl_data,options);
clear options;

% feedback QP controller with atlas
ins(1).system = 1;
ins(1).input = 2;
ins(2).system = 1;
ins(2).input = 3;
ins(3).system = 1;
ins(3).input = 4;
outs(1).system = 2;
outs(1).output = 1;
sys = mimoFeedback(qp,r_actual,[],[],ins,outs);
clear ins;

% feedback foot contact detector with QP/atlas
options.use_lcm=false;
options.contact_threshold = 0.002;
fc = FootContactBlock(r,ctrl_data,options);
ins(1).system = 2;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 3;
sys = mimoFeedback(fc,sys,[],[],ins,outs);
clear ins;

% feedback PD trajectory controller
options.use_ik = false;
pd = IKPDBlock(r,ctrl_data,options);
ins(1).system = 1;
ins(1).input = 1;
ins(2).system = 2;
ins(2).input = 2;
sys = mimoFeedback(pd,sys,[],[],ins,outs);
clear ins;

% set up a faked force/torque to one of the hands
rhand_idx = findLinkId(r,'r_hand');
% this could be some much better if we actually did the right math
ft = ConstantTrajectory([rhand_idx;0;0;0;-2;0;0]);
ft = ft.setOutputFrame(atlasFrames.ExternalForceTorque);
ins(1).system = 2;
ins(1).input = 1;
sys = mimoCascade(ft,sys,[],ins,outs);
clear ins;

qt = QTrajEvalBlock(r,ctrl_data);
sys = mimoFeedback(qt,sys,[],[],[],outs);

if visualize
  v = r.constructVisualizer;
  v.display_dt = 0.01;
  S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
  output_select(1).system=1;
  output_select(1).output=1;
  sys = mimoCascade(sys,v,[],[],output_select);
  warning(S);
end
x0 = xstar;
x0(3) = 1.0; % drop it a bit

traj = simulate(sys,[0 2],x0);
if visualize
  playback(v,traj,struct('slider',true));
end

xf = traj.eval(traj.tspan(2));

err = norm(xf(1:6)-xstar(1:6))
if err > 0.02
  error('drakeBalancing unit test failed: error is too large');
end

end
