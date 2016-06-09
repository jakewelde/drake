function runGSSphereMakeImages

% Make a world to collide against. It consists of a LEGO brick,
% like the reference data
r = GSSimpleSphere();
x0 = Point(getStateFrame(r));
use_margins = false;

refs_list = dir('spherereference');
% image_index = 72;
% ref_img = imread(['reference/',refs_list(2+image_index).name]);
% ref_img = double(ref_img) / max(max(max(double(ref_img))));

good_image_indices = [1, 65, 138, 290, 333, 334, 335, 336 ...
    417, 420, 421 ...
    ];

% Important globals
gs_pixelw = 200*4;
gs_pixelh = 150*4;
gs_width = .47*1.25;
gs_height = .37*1.25;
gs_thickness = 0.02;

% Hand-coded reference data:
transforms = {};
for i=1:max(good_image_indices)
    transforms{i} = [0,0,-10,0,0,0];
end
transforms{1} = [0,0,-10,0,0,0];
transforms{65} = [-0.155,0.232,-gs_thickness/2*.65,0,0,0];
transforms{138} = [-0.125,-0.12,-gs_thickness/2*1.1,0,0,0];
transforms{290} = [0.01,0.06,0.0,0,0,0];
transforms{333} = [0.1495,-0.1225,-gs_thickness/2,0,0,0];
transforms{334} = [0.147,-0.123,-gs_thickness/2*.7,0,0,0];
transforms{335} = [0.142,-0.122,-gs_thickness/2*.55,0,0,0];
transforms{336} = [0.137,-0.122,-gs_thickness/2*.37,0,0,0];

for i=1:max(good_image_indices)
    t = transforms{i};
    t(3) = t(3) + .00;
    transforms{i} = t;
end

for image_index=good_image_indices
    ref_img = imread(['spherereference/',refs_list(2+image_index).name]);
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

    xi = x0;

    [gs_manifold_x, gs_manifold_n]=makeGSRectangular(gs_pixelw, gs_pixelh, gs_width, gs_height, gs_thickness);

    gs_manifold_x = gs_manifold_x + repmat(approach_vec, 1, size(gs_manifold_x,2));

    % Position the box properly
    %xi.base_roll = xi.base_roll + pi/4;
    %xi.base_pitch = xi.base_pitch + .05;
    kinsol = r.doKinematics(double(xi));

    %(FOR NOW) Snap the sensor to the collision hull
    gs_manifold_x_unslided = gs_manifold_x;
    slide_amount = norm(approach_vec) - gs_thickness*2.5;

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

    compos_img = (.2*ref_img + .8*distances_rsz);

    imshow(compos_img);
    
    imwrite(ref_img,['cleanspherergb/ref' int2str(image_index) '.png']);
    imwrite(distances_rsz,['cleanspheredepth/ref' int2str(image_index) '.png']);

end

end