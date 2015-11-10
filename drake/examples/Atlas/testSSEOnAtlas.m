function findAtlasFixedPoint()

checkDependency('gurobi')

visualize = true;

% silence some warning
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

addpath('/home/gizatt/drc/software/drake/drake/solvers/dev');

options.floating = true;
options.dt = 0.002;
r = Atlas('/home/gizatt/drc/software/drake/drake/examples/Atlas/urdf/atlas_convex_hull.urdf',options);
r_bc = r;
r = r.removeCollisionGroupsExcept({'heel','toe', 'midfoot'});
%r = r.removeCollisionGroupsExcept({'heel','toe', 'midfoot', 'left_hand_robotiq', 'right_hand_robotiq'});
%r = r.removeCollisionGroupsExcept({'left_hand_robotiq', 'right_hand_robotiq'});
%r = r.removeCollisionGroupsExcept({'left_hand_robotiq', 'right_hand_robotiq', 'default'});
r = compile(r);

nq = getNumPositions(r);
nu = r.getNumInputs();
v = r.constructVisualizer;

joints = r.getStateFrame.getCoordinateNames();
joints = joints(1:nq);
l_arm_joints = find(strncmp(joints, 'l_arm_', 5));
r_arm_joints = find(strncmp(joints, 'r_arm_', 5));
n_arm_joints = numel(r_arm_joints);

l_leg_joints = find(strncmp(joints, 'l_leg_', 5));
r_leg_joints = find(strncmp(joints, 'r_leg_', 5));
n_leg_joints = numel(r_leg_joints);


% ytraj
l = load('walking_traj.mat');
x_gt = l.ytraj.eval(1.0);

% generate a point cloud
rad0 = pi/2;
rayAngles = linspace(-rad0, rad0, 400);
spindleAngles = linspace(-rad0, rad0, 400);
% position of the laser source
scanOffset = [-5; 0; 3];
% generate ray set
[theta_x, theta_z] = ndgrid(rayAngles, spindleAngles);
totalRange = 30;
raySet = totalRange*[cos(theta_x(:)), ...
    cos(theta_z(:)).*sin(theta_x(:)), ...
    sin(theta_z(:)).*sin(theta_x(:))]';
origin = repmat(scanOffset,1,numel(theta_x)); 
% and do raycast
q_gt = x_gt(1:nq);
kinsol_bc = r_bc.doKinematics(q_gt);
[distances, normals] = collisionRaycast(r_bc, kinsol_bc, origin, raySet+origin, 0);
pointsInWorldFrame = repmat(scanOffset,1,numel(distances)) + repmat(distances.', [3, 1]) .* raySet / totalRange;
pointsInWorldFrame = pointsInWorldFrame(:, distances >= 0 & (pointsInWorldFrame(3,:)>0.05).');

lcmgl = LCMGLClient('atlaspts_debug');
lcmgl.glColor3f(1,0,0);
lcmgl.points(pointsInWorldFrame(1,:), pointsInWorldFrame(2,:), pointsInWorldFrame(3,:));
lcmgl.switchBuffers();

scale_sequence = [0.00001;];

guess = struct();
guess.q = double(q_gt) + [0.1*(rand(6,1)-0.5); 0.25*(rand(numel(q_gt)-6, 1)-0.5)];

for i=1:numel(scale_sequence)
    scale = 1;
    options.compl_slack = scale*0.01;
    options.lincompl_slack = scale*0.001;
    options.jlcompl_slack = scale*0.01;
    options.scaling = scale_sequence(i);
    opt = ContactImplicitFixedPointUnconstrainedProgram(r.getManipulator(), [], options); 


    opt = opt.setSolverOptions('snopt','DerivativeOption', 0);
    opt = opt.setSolverOptions('snopt','VerifyLevel', 0);
    opt = opt.setSolverOptions('snopt','MajorOptimalityTolerance', 1E-5);
    opt = opt.setSolverOptions('snopt','SuperbasicsLimit', 10000);
    opt = opt.setSolverOptions('snopt','print','snopt.out');


    %opt = opt.addInputCost(QuadraticSumConstraint(0,0,0.1*eye(nu),zeros(nu,1)));
    % remove floating base fredom that we don't need
    %opt = opt.addConstraint(ConstantConstraint(zeros(3,1)), opt.q_inds([1 2 6]));
    %opt = opt.addConstraint(BoundingBoxConstraint(-0.1*ones(2,1), 0.1*ones(2,1)), opt.q_inds([4 5]));
    %opt = opt.addConstraint(BoundingBoxConstraint(0.6, 1.2), opt.q_inds([3]));
    %opt = opt.addCost(LinearConstraint(0, 0, -100), opt.q_inds([3]));

    % symmetry in the arms
    flip = [-1 -1 1 -1 1 -1 1];
    A = [eye(n_arm_joints) -diag(flip)];
    %opt = opt.addConstraint(LinearConstraint(zeros(n_arm_joints,1), zeros(n_arm_joints,1), A), [l_arm_joints;r_arm_joints]);

    % symmetry in the legs
    flip = [1 -1 1 1 1 -1];
    A = [eye(n_leg_joints) -diag(flip)];
    %opt = opt.addConstraint(LinearConstraint(zeros(n_leg_joints,1), zeros(n_leg_joints,1), A), [l_leg_joints;r_leg_joints]);

    % require contact force from feet
    %opt = opt.addConstraint(BoundingBoxConstraint(zeros(opt.nC,1)+1, zeros(opt.nC,1)+Inf), opt.l_inds(1:(opt.nD+1):end));
    %cp = 8;
    %opt = opt.addCost(FunctionHandleConstraint(0, 0, nq, @feet_on_ground_constraint_fun), {opt.q_inds});


    % INFORM THIS BY THE JOINT ENCODERS?
    % MAYBE THIS MAKES MORE SENSE IN LOW DIMENSION?
    % DUNNO...

    % with an SDF objective on those raycast points
    sdfCost = FunctionHandleConstraint(0,0,nq, @(q)staticStableEstimatorBox_SDFObjective(q,r_bc.getManipulator(),pointsInWorldFrame));
    %sdfCost.grad_method = 'numerical';
    opt = opt.addCost(sdfCost, opt.q_inds);

    v.draw(0, guess.q);
    [q, u, l, info, F] = opt.findFixedPoint(guess, v);
    guess.q = q;
    guess.u = u;
    guess.l = l;
    info
end

function [f,df] = feet_on_ground_constraint_fun(q)
    colopts = opt.options.active_collision_options;
    colopts.terrain_only = true;
    [phi,normal,~,~,~,idxA,idxB,~,n,D,dn,dD] = r.getManipulator().contactConstraints(q,opt.options.multiple_contacts,colopts);
    f = 100*phi.'*phi;
    df = 200*phi.'*n;
    df(:,1:6) = 0;
end

end
