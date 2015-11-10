options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.001;
%options.use_bullet = true;
%options.enable_fastqp = false;
options.use_mex = true;


w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile(getDrakePath,'systems','plants','test','FallingBrickContactPoints.urdf');
% and one that has an actual box for collision geometry:
urdf_bc = fullfile(getDrakePath,'systems','plants','test','FallingBrick.urdf'); 
p = RigidBodyManipulator(urdf, options);
p_bc = RigidBodyManipulator(urdf_bc, options);
ts = TimeSteppingRigidBodyManipulator(p, 0.001);
v = p.constructVisualizer();

nq = p.getNumPositions();

q_gt = [0 0 0.75 0 0 0].';
v.draw(0, q_gt);

% generate a point cloud
rad0 = pi/2;
rayAngles = linspace(-rad0, rad0, 200);
spindleAngles = linspace(-rad0, rad0, 200);
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
kinsol = p.doKinematics(q_gt);
kinsol_bc = p_bc.doKinematics(q_gt);
[distances, normals] = collisionRaycast(p_bc, kinsol_bc, origin, raySet+origin, 0);
pointsInWorldFrame = repmat(scanOffset,1,numel(distances)) + repmat(distances.', [3, 1]) .* raySet / totalRange;
pointsInWorldFrame = pointsInWorldFrame(:, distances >= 0 & (pointsInWorldFrame(3,:)>0.05).');

lcmgl = LCMGLClient('boxpts_debug');
lcmgl.glColor3f(1,0,0);
lcmgl.points(pointsInWorldFrame(1,:), pointsInWorldFrame(2,:), pointsInWorldFrame(3,:));
lcmgl.switchBuffers();

% Now set up our fixed point program
scale = 0.0;
options.compl_slack = scale*0.01;
options.lincompl_slack = scale*0.001;
options.jlcompl_slack = scale*0.01;
opt = ContactImplicitFixedPointUnconstrainedProgram(p, [], options); 
%opt = ContactImplicitFixedPointProgram(p, [], options); 
%opt = ContactImplicitLogForceFixedPointProgram(p, [], options); 
opt = opt.setSolverOptions('snopt','DerivativeOption', 0);
opt = opt.setSolverOptions('snopt','VerifyLevel', 0);
opt = opt.setSolverOptions('snopt','MajorOptimalityTolerance', 1E-5);
opt = opt.setSolverOptions('snopt','print','snopt.out');

% with an SDF objective on those raycast points
sdfCost = FunctionHandleConstraint(0,0,nq, @(q)staticStableEstimatorBox_SDFObjective(q,p_bc,pointsInWorldFrame));
%sdfCost.grad_method = 'numerical';
opt = opt.addCost(sdfCost, opt.q_inds);

guess = struct();
guess.q = q_gt + 0.05*(2*rand(6,1)-1);
[q, u, l, info, F] = opt.findFixedPoint(guess, v);
%q = opt.solve(guess.q);

fprintf('This is what the trajopt came up with as a physically reasonable\n');
fprintf('   fixed point near that cloud fit\n');
info
F