[terrainX, terrainY] = ndgrid([-5:0.1:5], [-5:0.1:5]);
% generate messy terrain with random + pink noise
freq = rand(size(terrainX));
Z = rand(size(terrainX))-0.5;
Z = fft2(Z);
alpha = 2;
K = 0.001;
dist = repmat([1:size(Z,1)].', 1, size(Z,2)) .* repmat(1:size(Z,2),size(Z,1), 1);
Z = Z ./ (K * dist.^alpha);
terrainZ = real(ifft2(Z));
% actually gonna use something easier that sstill breaks it
terrainZ = 0.1*terrainY.*terrainX;
%surf(terrainX, terrainY, terrainZ);

options.terrain = RigidBodyHeightMapTerrain(terrainX(:,1),terrainY(1,:),terrainZ);
options.floating = true;
options.dt = 0.001;
%options.use_bullet = true;
%options.enable_fastqp = false;
options.use_mex = true;


w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile(getDrakePath,'systems','plants','test','FallingBrick.urdf');
% and one that has an actual box for collision geometry:
urdf_bc = fullfile(getDrakePath,'systems','plants','test','FallingBrick.urdf'); 
p = RigidBodyManipulator(urdf, options);
p_bc = RigidBodyManipulator(urdf_bc, options);
ts = TimeSteppingRigidBodyManipulator(p, 0.001);
v = p.constructVisualizer();

nq = p.getNumPositions();

q_gt = [0 0 0.6 0 0 0].';
v.draw(0, q_gt);

% Now set up our fixed point program
scale = 0.0;
options.compl_slack = scale*0.01;
options.lincompl_slack = scale*0.001;
options.jlcompl_slack = scale*0.01;
opt = ContactImplicitFixedPointProgram(p, [], options); 
opt = opt.setSolverOptions('snopt','DerivativeOption', 0);
opt = opt.setSolverOptions('snopt','VerifyLevel', 0);
opt = opt.setSolverOptions('snopt','MajorOptimalityTolerance', 1E-5);
opt = opt.setSolverOptions('snopt','print','snopt.out');

opt = opt.addConstraint(BoundingBoxConstraint([-4;-4],[4;4]),opt.q_inds(1:2));

guess = struct();
guess.q = q_gt + 0.1*(2*rand(6,1)-1);
[q, u, l, info, F] = opt.findFixedPoint(guess, v);
%q = opt.solve(guess.q);

fprintf('This is what the trajopt came up with as a physically reasonable\n');
fprintf('   fixed point near that cloud fit\n');
info
F

if (info == 1)
    keyboard

    ts = TimeSteppingRigidBodyManipulator(urdf, 0.001, options);
    v_vis = v.setNumInputs(ts.getNumStates);
    v_vis = v_vis.setInputFrame(ts.getOutputFrame);
    output_select(1).system=1;
    output_select(1).output=1;
    sys = mimoCascade(ts,v_vis,[],[],output_select);
    'simulating to show that its a fixed pt'
    traj = sys.simulate([0,2.0], [q;zeros(6,1)]);
else
    'bad info, not simulating'
end