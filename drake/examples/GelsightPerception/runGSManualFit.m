function runGSManualFit

% Make a world to collide against. It consists of a LEGO brick,
% like the reference data
r = GSLegoEnd();
x0 = Point(getStateFrame(r));
use_margins = false;

refs_list = dir('reference');
image_index = 2+22;
ref_img = imread(['reference/',refs_list(image_index).name]);
ref_img = double(ref_img) / max(max(max(double(ref_img))));

% Move sensor down
dz = -10.0;
approach_vec = [0;0;dz];

% Important globals
gs_pixelw = 200;
gs_pixelh = 200;
gs_width = 1.75;
gs_height = 1.31;
gs_thickness = 0.2;

xi = x0;
%xi.base_roll = xi.base_roll + pi/4;

[gs_manifold_x, gs_manifold_n]=makeGSRectangular(gs_pixelw, gs_pixelh, gs_width, gs_height, gs_thickness);

gs_manifold_x = gs_manifold_x + repmat(approach_vec, 1, size(gs_manifold_x,2));

% Position the box properly
%xi.base_roll = xi.base_roll + pi/4;
%xi.base_pitch = xi.base_pitch + .05;
kinsol = r.doKinematics(double(xi));

%(FOR NOW) Snap the sensor to the collision hull
gs_manifold_x_unslided = gs_manifold_x;
slide_amount = collisionApproachGelSight(r, kinsol, gs_manifold_x_unslided, -approach_vec*2, use_margins);

slide_vector = (-approach_vec) * (slide_amount-(gs_thickness/2)) / norm(approach_vec);
gs_manifold_x = gs_manifold_x_unslided + repmat(slide_vector, 1, size(gs_manifold_x,2)); %.9945

% Collide! Make image from the rectangular GelSight sensor
distances = collisionGelSight(r, kinsol, gs_manifold_x, gs_manifold_n, use_margins);
distances = gs_thickness - distances;

% Reshape result and render image
distances = distances / gs_thickness;
distances = reshape(distances,[gs_pixelh, gs_pixelw]);

%Finally, show the output image
disp(size(ref_img))
imresize(distances, [size(ref_img,1),size(ref_img,2)]);
distances = repmat(distances,1,1,3);
disp(size(distances))
refw = size(ref_img,1);
refh = size(ref_img,2);
% YES , THIS IS OBVIOUSLY IN NEED OF REFACTORING
for x=0:(refw-1)
    for y=0:(refh-1)
        dx = floor(x * gs_pixelw / refw);
        dy = floor(y * gs_pixelh / refh);
        ref_img(x+1,y+1,:) = min(255, (ref_img(x+1,y+1,:) + 255*distances(dx+1,dy+1))/2);
    end
end

imshow(ref_img);

end