[terrainX, terrainY] = ndgrid([-5:0.1:5], [-5:0.1:5]);
% generate messy terrain with random + pink noise
freq = rand(size(terrainX));
terrain = 'tiltedblocks';
if strcmp(terrain, 'smoothsurf')
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
elseif strcmp(terrain, 'tiltedblocks')
    r = RigidBodyManipulator();
    r.setTerrain(RigidBodyFlatTerrain());
    
    rpy = [0;0;0];
    b = RigidBodyBox([5; 5; 0.1], [0;0;-0.05], zeros(3, 1));
    r = r.addGeometryToBody('world', b);
    for j = 1:100
      box_size = [0.39, 0.39, 0.146];
      rpy = [0;-15*pi/180;rand*5*pi];
      box_xyz = [rand(2, 1)*5-2.5; rand*0.5];
      b = RigidBodyBox(box_size, box_xyz + [0;0;-box_xyz(3)/2], rpy);
      r = r.addGeometryToBody('world', b);
    end
    r = r.compile();
    options.terrain = RigidBodyHeightMapTerrain.constructHeightMapFromRaycast(r,[],-2.5:.015:2.5,-2.5:.015:2.5,10);
else
    'invalid terrain'
    keyboard
end

options.floating = true;
options.dt = 0.001;
options.use_bullet = true;
%options.enable_fastqp = false;
options.use_mex = true;
options.multiple_contacts = true;


w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile(getDrakePath,'systems','plants','test','Capsule.urdf');
% and one that has an actual box for collision geometry:
urdf_bc = fullfile(getDrakePath,'systems','plants','test','Capsule.urdf'); 
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

opt = opt.addConstraint(BoundingBoxConstraint([-2.5;-2.5],[2.5;2.5]),opt.q_inds(1:2));

guess = struct();
guess.q = q_gt + 0.1*(2*rand(6,1)-1);
[q, u, l, info, F] = opt.findFixedPoint(guess, v);
%q = opt.solve(guess.q);

fprintf('This is what the trajopt came up with as a physically reasonable\n');
fprintf('   fixed point near that cloud fit\n');
info
F
% 
% if (info == 1)
%     keyboard
% 
%     ts = TimeSteppingRigidBodyManipulator(urdf, 0.001, options);
%     v_vis = v.setNumInputs(ts.getNumStates);
%     v_vis = v_vis.setInputFrame(ts.getOutputFrame);
%     output_select(1).system=1;
%     output_select(1).output=1;
%     sys = mimoCascade(ts,v_vis,[],[],output_select);
%     'simulating to show that its a fixed pt'
%     traj = sys.simulate([0,2.0], [q;zeros(6,1)]);
% else
%     'bad info, not simulating'
% end