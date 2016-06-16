function runGSMeasurementDemoRecon

    % Load measured images
    refs_list = dir('samplevideo');
    image_index = 1+187+0*112;
    
    % Important globals
    steps = 80;
    
    % Load and scale trained data
    tr = load('iamtrained.mat');
    %gs_bg = tr.a.gs_bg;
    gs_filter = tr.filter;
    scaling=.2;
    
    measured_img = double(imread(['samplevideo/',refs_list(4).name]));
    measured_img = measured_img / 255.0;
    measured_img = imresize(measured_img, scaling);
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
        measured_img = imresize(measured_img, scaling);%gs_pixelw/size(measured_img,2));
        
        ref_img_2D = measured_img - gs_bg;
        
        heightmap_recon_2D_qp = convertGStoHM(ref_img_2D,gs_filter,true);
        heightmap_recon_2D = convertGStoHM(ref_img_2D,gs_filter,false);
        
        %heightmap_recon_2D = (heightmap_recon_2D-min(min(min(heightmap_recon_2D)))) / ...
        %    (max(max(max(heightmap_recon_2D)))-min(min(min(heightmap_recon_2D))));
        
        %% Draw the reconstructed heightmap
        
        imshow(imresize(cat(2,ref_img_2D,repmat(heightmap_recon_2D_qp,[1 1 3]),repmat(heightmap_recon_2D,[1 1 3])),3.0));
        drawnow;
%         disp(heightmap_recon_2D(1,1));
    end

    disp('final time:');
    disp(cputime - timebefore);
    disp('final rate:');
    disp(counter / (cputime - timebefore));
end