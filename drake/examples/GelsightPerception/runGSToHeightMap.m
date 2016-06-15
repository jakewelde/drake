function runGSMeasurementDemoRecon

    % Load measured images
    refs_list = dir('samplevideo');
    image_index = 1+187+0*112;
    
    % Move sensor down
    dy = -10;
    approach_vec = [0;0;dy];

    % Important globals
    steps = 40;
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
% 
%     [gs_manifold_x_default, gs_manifold_n_default]=makeGSRectangular(gs_pixelw, gs_pixelh, gs_width, gs_height, gs_thickness);
% 
%     gs_manifold_x_default = gs_manifold_x_default + repmat(approach_vec, 1, size(gs_manifold_x_default,2));
% 
%     M = 150;

%     x_lo = [-gs_height/2,-gs_width/2,-2*gs_thickness,-pi/2,-pi/2,-pi/2];
%     x_hi = [gs_height/2,gs_width/2,0,pi/2,pi/2,pi/2];
%     x_sigma = (x_hi - x_lo)/6;
%     
%     parts = zeros(M,length(x_sigma));
%     for m=1:M
%         % Particles initialized to uniform distribution
%         parts(m,:) = rand(1,6).*(x_hi-x_lo) + x_lo;
%     end
%     parts_new = zeros(M,length(x_sigma)); % Allocate parts_new
%     weights_new = ones(M,1);
    
    % Load and scale trained data
    tr = load('iamtrained.mat');
    %gs_bg = tr.a.gs_bg;
    gs_filter = tr.filter;
    
    measured_img = double(imread(['samplevideo/',refs_list(4).name]));
    measured_img = measured_img / 255.0;
    measured_img = imresize(measured_img, gs_pixelw/size(measured_img,2));
    gs_bg = measured_img;
    
    timebefore = cputime;
    counter = 0;
    for t=1:steps
        counter = counter + 1;

        %% Grab next GelSight image
        image_index = image_index + 3;
        if (image_index > length(refs_list))
            disp('final time:');
            disp(cputime - timebefore);
            disp('final rate:');
            disp(counter / (cputime - timebefore))
            return;
        end
        measured_img = double(imread(['samplevideo/',refs_list(image_index).name]));
        measured_img = measured_img / 255.0;
        measured_img = imresize(measured_img, .2);%gs_pixelw/size(measured_img,2));
        
        ref_img_2D = measured_img - gs_bg;
        recon_2D = zeros(size(ref_img_2D));
        for color1=1:size(ref_img_2D,3)
            for color2=1:size(ref_img_2D,3)
                recon_2D(:,:,color2) = recon_2D(:,:,color2) + conv2(ref_img_2D(:,:,color1), gs_filter(:,:,color1,color2), 'same');
            end
        end

        hsize = floor(size(recon_2D,1)/10);
        sigma = floor(hsize/8);
        for color=1:3
            recon_2D(:,:,color) = conv2(recon_2D(:,:,color),fspecial('gaussian',hsize,sigma),'same');
        end
        
        % Use un-magnitude-corrected normals to reconstruct a heightmap
        normals_recon_2D = recon_2D;

%        assert(min(min(min(normals_recon_2D(:,:,3))))<=0, 'normal z comp 0 or less: %f',min(min(min(normals_recon_2D(:,:,3)))));
%        assert(min(min(min(normals_recon_2D(:,:,3))))<=0.1, 'normal z comp 0.1 or less: %f',min(min(min(normals_recon_2D(:,:,3)))));

        gradr_recon_2D = zeros(size(normals_recon_2D,1),size(normals_recon_2D,2));
        gradr_recon_2D(:,:) = -normals_recon_2D(:,:,1);% ./ normals_recon_2D(:,:,3);
        gradc_recon_2D = zeros(size(normals_recon_2D,1),size(normals_recon_2D,2));
        gradc_recon_2D(:,:) = -normals_recon_2D(:,:,2);% ./ normals_recon_2D(:,:,3);


        grads_recon_2D = 100*cat(3,gradr_recon_2D,gradc_recon_2D,zeros(size(gradr_recon_2D))); % + .0094;

        %% Build a QP to find a nice heightmap

        gradr_vec = cat(1,zeros(1,size(gradr_recon_2D,2)),gradr_recon_2D,zeros(1,size(gradr_recon_2D,2)));
        gradr_vec = gradr_vec(:);
        gradc_vec = cat(2,zeros(size(gradc_recon_2D,1),1),gradc_recon_2D,zeros(size(gradc_recon_2D,1),1));
        gradc_vec = gradc_vec(:);

        %Sparse matrices encoding the row- and col-wise gradient operations
        Dr = convmtx2([1;0;-1],[size(grads_recon_2D,1),size(grads_recon_2D,2)]);
        Dc = convmtx2([1,0,-1],[size(grads_recon_2D,1),size(grads_recon_2D,2)]);

        border_left_2D = cat(2,ones(size(grads_recon_2D,1),1),zeros(size(grads_recon_2D,1),size(grads_recon_2D,2)-1))==1;
        border_top_2D = cat(1,ones(1,size(grads_recon_2D,2)),zeros(size(grads_recon_2D,1)-1,size(grads_recon_2D,2)))==1;
        border_full_2D = (0*grads_recon_2D(:,:,1))==1; %initially full of false
        border_full_2D = border_full_2D | border_left_2D; %set one border
        border_full_2D = border_full_2D | border_top_2D; %set one border
        border_full_2D = border_full_2D | imrotate(border_left_2D,180); %set one border
        border_full_2D = border_full_2D | imrotate(border_top_2D,180); %set one border
        border_full_2D = double(border_full_2D);
        Q = spdiags(border_full_2D(:),[0],numel(border_full_2D),numel(border_full_2D));

        H = Dr'*Dr + Dc'*Dc + Q'*Q;
        f = -2*gradr_vec'*Dr + -2*gradc_vec'*Dc;
        
        x0 = 0*gradr_recon_2D(:);
        x0 = x0 + .001;

        opts = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','none');
%         opts = optimoptions('quadprog',...
%         'Algorithm','interior-point-convex','Display','iter');
        heightmap_recon = quadprog(H,f,[],[],[],[],[],[],[],opts);

        heightmap_recon_2D = reshape(heightmap_recon,size(gradr_recon_2D));
        heightmap_recon_2D = max(0,heightmap_recon_2D);
        heightmap_recon_2D = 300*heightmap_recon_2D;
        %heightmap_recon_2D = (heightmap_recon_2D-min(min(min(heightmap_recon_2D)))) / ...
        %    (max(max(max(heightmap_recon_2D)))-min(min(min(heightmap_recon_2D))));
        %% Draw the reconstructed heightmap
        
        disp(heightmap_recon_2D(1,1));
%         imshow(imresize(cat(2,ref_img_2D,repmat(heightmap_recon_2D>.5,[1 1 3])),3.0));
%         drawnow;
        
%         % Amplify those particles which have higher weight
%         sampled_indices = randsample(M,M,true,weights_new);
%         parts = parts_new(sampled_indices,:);
 
%         figure
%         scatter(parts(:,1),parts(:,2));
    end

%     function distances=make_image(x)
%         % Position the box properly
%         xi = x0;
%         xi.base_x = x(1);
%         xi.base_y = x(2);
%         %xi.base_z = x(3); %% See below
%         xi.base_roll = x(4);
%         xi.base_pitch = x(5);
%         xi.base_yaw = x(6);
%         kinsol = r.doKinematics(double(xi));
% 
%         % Snap the sensor to the collision hull
%         gs_manifold_x_unslided = gs_manifold_x_default;
%         slide_amount = collisionApproachGelSight(r, kinsol, gs_manifold_x_unslided, -approach_vec*2, use_margins);
% 
%         slide_vector = (-approach_vec) * (slide_amount) / norm(approach_vec)  +  x(3);
%         gs_manifold_x = gs_manifold_x_unslided + repmat(slide_vector, 1, size(gs_manifold_x_unslided,2));
%         gs_manifold_n = gs_manifold_n_default;
% 
%         % Collide! Make image from the rectangular GelSight sensor
%         distances = collisionGelSight(r, kinsol, gs_manifold_x, gs_manifold_n, use_margins);
%         distances = gs_thickness - distances;
% 
%         % Reshape result and render image
%         distances = distances / gs_thickness;
%         distances = reshape(distances,[gs_pixelh, gs_pixelw]);
%         distances = repmat(distances,[1,1,3]);
%     end
    disp('final time:');
    disp(cputime - timebefore);
    disp('final rate:');
    disp(counter / (cputime - timebefore));
end