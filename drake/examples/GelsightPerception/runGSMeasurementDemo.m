function runGSMeasurementDemo

    % Load measured images
    refs_list = dir('samplevideo');
    image_index = 1+187+0*112;
    
    % Make a world to collide against
    r = GSLegoEndSpecial();% GSSimpleBox();
    % v = r.constructVisualizer();
    x0 = Point(getStateFrame(r));
    use_margins = false;

    % Move sensor down
    dy = -10;
    approach_vec = [0;0;dy];

    % Important globals
    steps = 4000;
    scaling = 1.0; % image scaling, for performance
    
    gs_pixelw = 128;
    %gs_width = 1.0;
    ar = 3/4;
    gs_pixelh = floor(ar*gs_pixelw);
    gs_width = 1.88;
    gs_height = 1.48;
    gs_thickness = 0.2;

    % Set up fast canvas
    [gs_manifold_x_default, gs_manifold_n_default]=makeGSRectangular(gs_pixelw, gs_pixelh, gs_width, gs_height, gs_thickness);
    gs_manifold_x_default = gs_manifold_x_default + repmat(approach_vec, 1, size(gs_manifold_x_default,2));

    M = 150;

    x_lo = [-gs_height/2,-gs_width/2,-gs_thickness,-pi/2,-pi/2,-pi/2];
    x_hi = [gs_height/2,gs_width/2,0,pi/2,pi/2,pi/2];
    x_sigma = (x_hi - x_lo)/6;
    
    parts = zeros(M,length(x_sigma));
    for m=1:M
        % Particles initialized to uniform distribution
        parts(m,:) = rand(1,6).*(x_hi-x_lo) + x_lo;
    end
    parts_new = zeros(M,length(x_sigma)); % Allocate parts_new
    weights_new = ones(M,1);
    
    % Load and scale trained data
    tr = load('iamtrained.mat');
    gs_filter = tr.filter;
    
    measured_img = double(imread(['samplevideo/',refs_list(4).name]));
    measured_img = measured_img / 255.0;
    measured_img = imresize(measured_img, gs_pixelw/size(measured_img,2));
    gs_bg = measured_img;
    
    one_true_x = parts(1,:);
    
    for t=1:steps

        % Grab next GelSight image
        image_index = image_index + 3;
        assert(image_index <= length(refs_list));
        measured_img = double(imread(['samplevideo/',refs_list(image_index).name]));
        measured_img = measured_img / 255.0;
        measured_img = imresize(measured_img, scaling * gs_pixelw/size(measured_img,2));
        measured_img = measured_img - gs_bg;
        
        heightmap_recon_2D = convertGStoHM(measured_img,gs_filter);
        heightmap_recon_2D = sqrt(max(0,heightmap_recon_2D));
        
        %DEBUGGGG
        heightmap_recon_2D = make_image(one_true_x);
        heightmap_recon_2D = imresize(heightmap_recon_2D,scaling);
        
        disp(size(heightmap_recon_2D));
        
        depth_imgs = cell(M,1);
        
        % Mutate current batch of particles
        for m=1:M
            % Gaussian-sample, corrected to be in-bounds
            parts_new(m,:) = max(x_lo, min(x_hi, normrnd(parts(m,:),x_sigma)));
            
            % Derive weight through image
            depth_img_2D = make_image(parts_new(m,:));
            depth_img_2D = imresize(depth_img_2D, scaling);
            
            depth_imgs{m,1} = depth_img_2D;

            weights_new(m) = probMatchGelSightDepth(depth_img_2D, heightmap_recon_2D);
        end
        
        B = zeros(M,M);
        for i=1:M
            for j=1:M
                B(i,j) = probMatchGelSightDepth(depth_imgs{i,1},depth_imgs{j,1});
            end
        end
        disp(B(1:5,1:5));
                
        %Some nice pictures
        [~,best_ind] = max(weights_new);
        mean_x = parts_new(best_ind,:);
        %mean_x = weights_new' * parts_new;
        mean_depth_img = make_image(mean_x);
        imshow(imresize(cat(2,repmat(mean_depth_img,[1 1 3]),repmat(heightmap_recon_2D,[1 1 3]),measured_img),3.0));
        %surf(imrotate(heightmap_recon_2D,180));
        %axis([0 size(heightmap_recon_2D,2) 0 size(heightmap_recon_2D,1) 0 1]);
        drawnow;
        
        % Amplify those particles which have higher weight
        sampled_indices = randsample(M,M,true,weights_new);
        parts = parts_new(sampled_indices,:);
        
%         figure
%         scatter(parts(:,1),parts(:,2));
    end

    function distances=make_image(x)
        % Position the box properly
        xi = x0;
        xi.base_x = x(1);
        xi.base_y = x(2);
        %xi.base_z = x(3); %% See below
        xi.base_roll = x(4);
        xi.base_pitch = x(5);
        xi.base_yaw = x(6);
        kinsol = r.doKinematics(double(xi));

        % Snap the sensor to the collision hull
        gs_manifold_x_unslided = gs_manifold_x_default;
        slide_amount = collisionApproachGelSight(r, kinsol, gs_manifold_x_unslided, -approach_vec*2, use_margins);

        slide_vector = (-approach_vec) * (slide_amount) / norm(approach_vec)  +  x(3);
        gs_manifold_x = gs_manifold_x_unslided + repmat(slide_vector, 1, size(gs_manifold_x_unslided,2));
        gs_manifold_n = gs_manifold_n_default;

        % Collide! Make image from the rectangular GelSight sensor
        distances = collisionGelSight(r, kinsol, gs_manifold_x, gs_manifold_n, use_margins);
        distances = gs_thickness - distances;

        % Reshape result and render image
        distances = distances / gs_thickness;
        distances = reshape(distances,[gs_pixelh, gs_pixelw]);
    end

end