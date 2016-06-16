function runGSToHeightMap(render_display)

    if nargin < length({'render_display'})
        render_display = false;
    end
    
    if render_display
        figure;
    end
    
    % Load trained filter
    tr = load('iamtrained.mat');
    gs_filter = tr.filter;
    
    gs_bg = 0; % Will pull first gelsight image as background
    gs_bg_initialized = false;
    
    lc = lcm.lcm.LCM.getSingleton();
    lcmonitor_image_raw = drake.util.MessageMonitor(bot_core.image_t,'utime');
    lc.subscribe('GELSIGHT_RAW',lcmonitor_image_raw);
    
    last_few_times = repmat(cputime, [5 1]);   
    last_few_times_marker = 1
    while true

        %% Grab next GelSight image
        data_image_raw = lcmonitor_image_raw.getMessage();
        if (~isempty(data_image_raw))
            data_image_t = bot_core.image_t(data_image_raw);
            data = data_image_t.data;
            
            width = data_image_t.width;
            height = data_image_t.height;
            
            measured_img = reshape(double(data),[3 width height]);
            measured_img = shiftdim(measured_img,1);
            measured_img = permute(measured_img,[2 1 3]);
            measured_img = measured_img / 255.0;
            
            if gs_bg_initialized ~= true
                gs_bg = measured_img;
                gs_bg_initialized=true;
                continue
            end
            
            ref_img_2D = measured_img - gs_bg;

            heightmap_recon_2D = convertGStoHM(ref_img_2D,gs_filter);        

            heightmap_recon_2D = heightmap_recon_2D';
            
            %% Publish
            lcm_depth_msg = bot_core.image_t();
            lcm_depth_msg.utime = cputime*1000*1000;
            lcm_depth_msg.width = size(heightmap_recon_2D,1);
            lcm_depth_msg.height = size(heightmap_recon_2D,2);
            lcm_depth_msg.row_stride = size(heightmap_recon_2D,1);
            lcm_depth_msg.pixelformat = 1497715271;
                % bot_core.image_t.BOT_CORE_IMAGE_T_PIXEL_FORMAT_GRAY;
            lcm_depth_msg.size = numel(heightmap_recon_2D);
            lcm_depth_msg.data = java.util.Arrays.copyOf(uint8(255*(heightmap_recon_2D(:))),length(heightmap_recon_2D(:)));
            lcm_depth_msg.nmetadata = 0;
            lc.publish('GELSIGHT_DEPTH', lcm_depth_msg);

            %% Draw the reconstructed heightmap

            %disp(heightmap_recon_2D(1,1));
            if render_display
                imshow(cat(2,ref_img_2D,repmat(heightmap_recon_2D',[1 1 3])));
                drawnow;
            end
            
            disp('Rate: ')
            disp(length(last_few_times)/(cputime-last_few_times(last_few_times_marker)));
            last_few_times(last_few_times_marker) = cputime;
            last_few_times_marker = mod(last_few_times_marker,length(last_few_times)) + 1;
            
        else
            % If image empty, wait a moment before continuing to avoid
            % thrashing
            pause(1/40);
        end
    end
end