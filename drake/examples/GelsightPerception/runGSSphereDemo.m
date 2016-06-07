function runGSSphereDemo

% Make a world to collide against
r = GSSimpleSphere();
x0 = Point(getStateFrame(r));
use_margins = false;
sphere_radius = .05;

% Move sensor down
dy = -10;
approach_vec = [0;0;dy];

% Important globals
steps = 100;
gs_pixelw = 200;
gs_width = 1.0;
gs_thickness = 0.005;

% Set up fast canvas
figure('doublebuffer','on');
colormap('Gray');
set(gca,'drawmode','fast');
set(gca,'units','pixels');
set(gca,'xlim',[0 gs_pixelw]);
set(gca,'ylim',[0 gs_pixelw]);
axis off;
imh = image('cdata', zeros(gs_pixelw, gs_pixelw));
set(imh,'erasemode','none');

x0.base_z = gs_thickness/2;
xi = x0;
for i=1:steps

    [gs_manifold_x, gs_manifold_n]=makeGSSquare(gs_pixelw, gs_width, gs_thickness);

    gs_manifold_x = gs_manifold_x + repmat(approach_vec, 1, size(gs_manifold_x,2));

    % Position the ball properly (slide it around, don't let it go too far)
    xi.base_xdot = xi.base_xdot + .01*(rand()-.5);
    xi.base_ydot = xi.base_ydot + .01*(rand()-.5);
    xi.base_zdot = xi.base_zdot + .001*(rand()-.5);
    if mod(i,ceil(steps/10)) == 0
        xi.base_xdot = 0;
        xi.base_ydot = 0;
        xi.base_zdot = 0;
    end
    xi.base_x = xi.base_x + xi.base_xdot;
    xi.base_y = xi.base_y + xi.base_ydot;
    xi.base_z = xi.base_z + xi.base_zdot;
    xi.base_x = max(-.4*gs_width, min(.4*gs_width, xi.base_x));
    xi.base_y = max(-.4*gs_width, min(.4*gs_width, xi.base_y));
    xi.base_z = max(-gs_thickness*.9, min(gs_thickness*.9, xi.base_z));
    kinsol = r.doKinematics(double(xi));
    
    %(FOR NOW) Snap the sensor to the collision hull
    gs_manifold_x_unslided = gs_manifold_x;
    slide_amount = norm(approach_vec) - sphere_radius;
    %slide_amount = collisionApproachGelSight(r, kinsol, gs_manifold_x_unslided, -approach_vec*2, use_margins);
    
    slide_vector = (-approach_vec) * (slide_amount) / norm(approach_vec);
    gs_manifold_x = gs_manifold_x_unslided + repmat(slide_vector, 1, size(gs_manifold_x,2)); %.9945

    % Collide! Make image from the rectangular GelSight sensor
    distances = collisionGelSight(r, kinsol, gs_manifold_x, gs_manifold_n, use_margins);
    distances = gs_thickness - distances;

    % Reshape result and render image
    distances = distances / gs_thickness;
    distances = reshape(distances,[gs_pixelw, gs_pixelw]);
    distances = uint8(distances * 255);
    % image(repmat(distances,1,1,3));
    set(imh,'cdata',distances);
    drawnow;
end

close all

end