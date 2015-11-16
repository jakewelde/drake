checkDependency('gurobi');
checkDependency('lcmgl');

robot_options = struct();
walking_options = struct();

robot_options = applyDefaults(robot_options, struct('use_bullet', true,...
                                                    'terrain', RigidBodyFlatTerrain,...
                                                    'floating', true,...
                                                    'ignore_self_collisions', true,...
                                                    'ignore_friction', true,...
                                                    'enable_fastqp', false,...
                                                    'dt', 0.001));
% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

% construct robot model
r_collision = Atlas('/home/gizatt/drc/software/drake/drake/examples/Atlas/urdf/atlas_convex_hull.urdf',robot_options);
r = Atlas(fullfile(getDrakePath,'examples','Atlas','urdf','atlas_minimal_contact.urdf'),robot_options);
%r = r.removeCollisionGroupsExcept({'heel','toe'});
r = compile(r);
% Instrument it
pitch = 0.5;
yaw = 0.5;
rows = 100;
cols = 100;
scanOffset = [1.5;0;0.5;0;0;pi];
rgbd_frame = RigidBodyFrame(1,scanOffset(1:3),scanOffset(4:6),'rgbdframe');
r_collision = r_collision.addFrame(rgbd_frame);
rgbd_sensor = RigidBodyDepthSensor('rgbd', r_collision.findFrameId('rgbdframe'), -pitch, pitch, rows, -yaw, yaw, cols, 10.0);
%rgbd_sensor = rgbd_sensor.enableLCMGL();
r_collision = r_collision.addSensor(FullStateFeedbackSensor());
r_collision = r_collision.addSensor(rgbd_sensor);
r_collision = r_collision.compile();


walking_options = applyDefaults(walking_options, struct('initial_pose', [],...
                                                        'navgoal', [1.5;0;0;0;0;0],...
                                                        'max_num_steps', 6,...
                                                        'rms_com_tolerance', 0.0051));
walking_options = applyDefaults(walking_options, r.default_footstep_params);
walking_options = applyDefaults(walking_options, r.default_walking_params);

% set initial state to fixed point
load(fullfile(getDrakePath,'examples','Atlas','data','atlas_fp.mat'));
if ~isempty(walking_options.initial_pose), xstar(1:6) = walking_options.initial_pose; end
xstar = r.resolveConstraints(xstar);
r = r.setInitialState(xstar);

v = r.constructVisualizer;
v.display_dt = 0.01;

nq = getNumPositions(r);

x0 = xstar;

% Find the initial positions of the feet
R=rotz(walking_options.navgoal(6));

rfoot_navgoal = walking_options.navgoal;
lfoot_navgoal = walking_options.navgoal;

rfoot_navgoal(1:3) = rfoot_navgoal(1:3) + R*[0;-0.13;0];
lfoot_navgoal(1:3) = lfoot_navgoal(1:3) + R*[0;0.13;0];

% Plan footsteps to the goal
goal_pos = struct('right', rfoot_navgoal, 'left', lfoot_navgoal);
footstep_plan = r.planFootsteps(x0(1:nq), goal_pos, [], struct('step_params', walking_options));
for j = 1:length(footstep_plan.footsteps)
  footstep_plan.footsteps(j).walking_params = walking_options;
end

% Generate a dynamic walking plan
walking_plan_data = r.planWalkingZMP(x0(1:r.getNumPositions()), footstep_plan);

options = struct();
options = applyDefaults(options, struct('gui_control_interface', false));

typecheck(r, 'Atlas');

if ~isfield(options, 'v') || isempty(options.v)
  v = r.constructVisualizer;
  v.display_dt = 0.05;
else
  v = options.v;
end

x0 = r.getInitialState();

% Build our controller and plan eval objects
control = atlasControllers.InstantaneousQPController(r, []);
planeval = atlasControllers.AtlasPlanEval(r, walking_plan_data);
plancontroller = atlasControllers.AtlasPlanEvalAndControlSystem(r, control, planeval);
sys = feedback(r_collision, plancontroller);

% Construct our QP estimator
qplc_est = QPLCContactImplicitEstimator(r_collision, x0, []);
clear outs;
for i=1:sys.getOutputFrame.getNumFrames
   outs(i).system = 1;
   outs(i).output = i;
end
sys = mimoCascade(sys, qplc_est, [], [], outs);

% Add a visualizer
output_select(1).system=1;
output_select(1).output=1;
sys = mimoCascade(sys,v,[],[],output_select);

% Simulate and draw the result
T = walking_plan_data.duration + 1;
x0_full = [x0; x0];
ytraj = simulate(sys, [0, T], x0_full, options);
