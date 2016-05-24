function runGSSim

close all

% Make a standard rectangular GS sensor, xyz-centered on the origin
w=200;
width=1.0;
gs_thickness = 0.05;
row = width * ( floor(-w/2)+1:floor(w/2) ) / w; % a row [-w, ..., w]
grid = cat(1,repmat(row,1,w),kron(row,ones(1,w))); % a 2-by-w*w mat, row x row:
 % [[-w, ..., w,   -w, ..., w,  .......  , -w, ..., w]
 %  [-w, ...,-w, -w+1,...,-w+1, .......  ,  w, ..., w]]
gs_manifold_x = cat(1,grid,zeros(1,w*w)); % 3-by-w*w, col-vecs of positions
gs_manifold_n = repmat([0;0;gs_thickness], 1, w*w); % 3-by-w*w, col-vecs of norms

% Make a world to collide against
r = GSSimpleBox();
% v = r.constructVisualizer();
kinsol = doKinematics(r, zeros(6,1));
use_margins = false;

%(FOR NOW) Move sensor down, then rotate it & normals about x-axis by tau/8
dy = -10; theta = pi/4;
R = eye(3);
%R = [1,0,0;0,cos(theta),-sin(theta);0,sin(theta),cos(theta)]'*R; %x-rot
R = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]'*R; %z-rot
R = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)]'*R; % y-rot
approach_vec_unrot = [0;0;dy];
approach_vec = R*approach_vec_unrot;

gs_manifold_x = gs_manifold_x + repmat(approach_vec_unrot, 1, w*w);
gs_manifold_x = R*gs_manifold_x;
gs_manifold_n = R*gs_manifold_n;

%(FOR NOW) Snap the sensor to the collision hull
gs_manifold_x_unslided = gs_manifold_x;
slide_amount = collisionApproachGelSight(r, kinsol, gs_manifold_x_unslided, -approach_vec*2, use_margins);
steps = 50;

for i=1:steps
    slide_vector = (-approach_vec) * (slide_amount-gs_thickness+(i*gs_thickness/steps)) / norm(approach_vec);
    gs_manifold_x = gs_manifold_x_unslided + repmat(slide_vector, 1, w*w); %.9945

    % Collide! Make image from the rectangular GelSight sensor
    distances = collisionGelSight(r, kinsol, gs_manifold_x, gs_manifold_n, use_margins);
    distances = gs_thickness - distances;

    % Reshape result and render image
    %%% TODO distances = (distances - min(distances)) / (max(distances) - min(distances));
    distances = distances / gs_thickness;
    distances = reshape(distances,[w,w]);
    image(repmat(distances,1,1,3));
    drawnow;
    hold on;
end

    function new_point=snap_point(rbm, kinsol, old_point)
        % old_point should be 3-by-1 point in global frame
        phi,n,x,body_x,body_idx = collisionDetectFromPoints(rbm,kinsol,old_point,false);
        new_point = x;
    end

end