function runGSMeasurementDemo

    % Load measured images
    refs_list = dir('samplevideo');
    image_index = 1+187+0*112;
    
    % Make a world to collide against
    r = GSLegoEndSpecial();% GSSimpleBox();
    v = r.constructVisualizer();
    x0 = Point(getStateFrame(r));
    use_margins = false;

    % Move sensor down
    dy = -10;
    approach_vec = [0;0;dy];

    % Important globals
    steps = 4000;
    gs_pixelw = 128;
    %gs_width = 1.0;
    ar = 3/4;
    gs_pixelh = floor(ar*gs_pixelw);
    %gs_height = ar*gs_width;
    %gs_thickness = 0.05;
    
    gs_width = 1.88;
    gs_height = 1.48;
    gs_thickness = 0.2;

    % Set up fast canvas
%     figure('doublebuffer','on');
%     colormap('Gray');
%     set(gca,'drawmode','fast');
%     set(gca,'units','pixels');
%     set(gca,'xlim',[0 gs_pixelw]);
%     set(gca,'ylim',[0 gs_pixelh]);
%     axis off;
%     imh = image('cdata', zeros(gs_pixelw, gs_pixelh));
%     set(imh,'erasemode','none');

    [gs_manifold_x_default, gs_manifold_n_default]=makeGSRectangular(gs_pixelw, gs_pixelh, gs_width, gs_height, gs_thickness);

    gs_manifold_x_default = gs_manifold_x_default + repmat(approach_vec, 1, size(gs_manifold_x_default,2));

    M = 150;

    x_lo = [-gs_height/2,-gs_width/2,-2*gs_thickness,-pi/2,-pi/2,-pi/2];
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
    tr = load('trained.mat');
    gs_bg = tr.a.gs_bg;
    gs_filter = tr.a.gs_filter;
    
    measured_img = double(imread(['samplevideo/',refs_list(4).name]));
    measured_img = measured_img / 255.0;
    measured_img = imresize(measured_img, gs_pixelw/size(measured_img,2));
    gs_bg = measured_img;
    
    for t=1:steps

        % Grab next GelSight image
        image_index = image_index + 3;
        assert(image_index <= length(refs_list));
        measured_img = double(imread(['samplevideo/',refs_list(image_index).name]));
        measured_img = measured_img / 255.0;
        measured_img = imresize(measured_img, gs_pixelw/size(measured_img,2));
        
        % Mutate current batch of particles
        for m=1:M
            % Gaussian-sample, corrected to be in-bounds
            parts_new(m,:) = max(x_lo, min(x_hi, normrnd(parts(m,:),x_sigma)));
            
            % Derive weight through image
            depth_img_2D = make_image(parts_new(m,:));
            
%             if m==1
%                 imshow(depth_img_2D);
%                 drawnow;
%                 figure;
%                 imshow(measured_img);
%                 drawnow;
%             end

            weights_new(m) = probMatchGelSightNaive(depth_img_2D, measured_img, gs_bg);
        end
                
        %Some nice pictures
        %[~,best_ind] = max(weights_new);
        %mean_x = parts_new(best_ind,:);
        mean_x = weights_new' * parts_new;
        mean_depth_img = make_image(mean_x);
%         mean_color_img = zeros(size(mean_depth_img,1), size(mean_depth_img,2), 1);
%         for color=1:3
%             mean_color_img(:,:,color) = conv2(mean_depth_img(:,:,color), gs_filter(:,:,color), 'same');
%         end
%         mean_color_img = mean_color_img + gs_bg;
        gradr_img_2D = conv2(mean_depth_img(:,:,1),[1,-1]','same');
        gradc_img_2D = conv2(mean_depth_img(:,:,1),[1,-1],'same');
        normals_img_2D = zeros(size(mean_depth_img));
        normals_img_2D(:,:,1) = -gradr_img_2D;
        normals_img_2D(:,:,2) = -gradc_img_2D;
        normals_img_2D(:,:,3) = ones(size(mean_depth_img(:,:,3)));
        norms_img_2D = zeros(size(normals_img_2D(:,:,1)));
        for direction=1:3
            norms_img_2D = norms_img_2D + normals_img_2D(:,:,direction) .* ...
                normals_img_2D(:,:,direction);
        end
        for direction=1:3
            normals_img_2D(:,:,direction) = normals_img_2D(:,:,direction) ./ ...
                norms_img_2D;
        end
        mags_img_2D = repmat(sin(2*acos(normals_img_2D(:,:,3))),[1,1,3]);
        mean_color_img = mags_img_2D;
        imshow(imresize(cat(2,mean_depth_img,mean_color_img,measured_img-gs_bg),3.0));
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
        distances = repmat(distances,[1,1,3]);
    end

end