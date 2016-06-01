function runGSManualFit

% Make a world to collide against. It consists of a LEGO brick,
% like the reference data
r = GSLegoEnd();
x0 = Point(getStateFrame(r));
use_margins = false;

refs_list = dir('reference');
% image_index = 72;
% ref_img = imread(['reference/',refs_list(2+image_index).name]);
% ref_img = double(ref_img) / max(max(max(double(ref_img))));

good_image_indices = [1,2,3,4, ...
    5,6,7,8, ...
    21,22,23,24, ...
    69,70,71,72, ...
    ]; %#ok<NASGU>

% Hand-coded reference data:
transforms = { ...
    [0.62,-0.836,0,0,0,.032], [0.62,-0.836,0,0,0,.032], [0.62,-0.836,0,0,0,.032], [0.62,-0.836,0,0,0,.032], ...
    [0.685,-.2,0,0,0,.04], [0.685,-.2,0,0,0,.04], [0.685,-.2,0,0,0,.04], [0.685,-.2,0,0,0,.04], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ... %20
    [0.1,-.175,0,0,0,.038], [0.1,-.175,0,0,0,.038], [0.1,-.178,0,0,0,.04], [0.1,-.178,0,0,0,.04], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ... %40
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ... %60
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0.36,-0.3,0,.9*pi/4,pi/10,.355*pi], [0.36,-0.3,0,.9*pi/4,pi/10,.355*pi], ...
    [0.36,-0.3,0,.9*pi/4,pi/10,.365*pi], [0.36,-0.3,0,.9*pi/4,pi/10,.365*pi], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ... %80
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ... %100
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], ...
    [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]};

for image_index=good_image_indices
    ref_img = imread(['reference/',refs_list(2+image_index).name]);
    ref_img = double(ref_img) / max(max(max(double(ref_img))));


    % Apply reference transform to brick
    transform = transforms{image_index};
    x0.base_x = transform(1);
    x0.base_y = transform(2);
    x0.base_z = transform(3);
    x0.base_roll = transform(4);
    x0.base_pitch = transform(5);
    x0.base_yaw = transform(6);

    % Move sensor down
    dz = -10.0;
    approach_vec = [0;0;dz];

    % Important globals
    gs_pixelw = 200*4;
    gs_pixelh = 150*4;
    gs_width = 1.88;
    gs_height = 1.48;
    gs_thickness = 0.2;

    xi = x0;

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
    % imresize(distances, [size(ref_img,1),size(ref_img,2)]);
    distances = repmat(distances,1,1,3);
    distances_rsz = imresize(distances, size(ref_img,1)/size(distances,1));

    disp(max(max(max(distances_rsz))));
    disp(max(max(max(ref_img))));

    compos_img = (.4*ref_img + .6*distances_rsz);

    imshow(compos_img);
    
    imwrite(ref_img,['cleanrgb/ref' int2str(image_index) '.png']);
    imwrite(distances_rsz,['cleandepth/ref' int2str(image_index) '.png']);

end

end