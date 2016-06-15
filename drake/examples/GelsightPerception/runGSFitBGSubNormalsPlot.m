function [filter, bg] = runGSFitBGSubNormalsPlot(full_model_lin)
% runGSFitBGSubGradient uses depthmaps and corresponding rgb
% images to train a linear convolution filter to map from heightmap to
% GelSight image.
%
% The linear model is as follows:
%
%   red_pixel = red_row_filter' * grad_row_pixels + ...
%                   red_col_filter' * grad_col_pixels;
%
% Where grad_row/col_pixels are a patch of row-/column-wise gradients
% derived from a depth image. There are corresponding, independent models
% for green and blue. By using raw gradients instead of the original depth
% map, we add some smoothness under different scalings of the filter (ie.
% one can first use the proper gradient operator on a high-res image, THEN
% apply the scaled filter to said image, rather than scaling both the
% gradient operator and the filter as a depth filter would require).
%
% @retval filter the flipped 2D filter image, m-by-n-by-3-by-2, with
% separate color channels for r,g, and b stored in its third dimension and
% separate directional gradient channels (1=row, 2=col) stored in its
% fourth dimension. Can be applied to a depth image to produce one
% color=1/2/3 of an rgb image by means of the MATLAB command conv2, eg:
%
%    grad_row_image(:,:,color)=conv2(depth_image(:,:,color),[1,-1]','same');
%    grad_col_image(:,:,color)=conv2(depth_image(:,:,color),[1,-1],'same');
%    gs_image(:,:,color) = ...
%      conv2(grad_row_image(:,:,color), filter(:,:,color,1), 'same') + ...
%      conv2(grad_col_image, filter(:,:,color,2), 'same');
%
% @retval bg the background that was subtracted from the images. Helpful
% for visualizing the derived GelSight image from a depth map, also shows
% what size of image the filter was derived from.

if nargin < length({'desired_img_width'})
   desired_img_width = 1280; 
end

refs_list = dir('reference');

% cleanrgbfolder = 'cleanrgb';
% cleandepthfolder = 'cleandepth';
% background_image_index = 1;
% good_image_indices = [1, ...
%     5, ...
%     21, ...
%     69 ...
%     ];

cleanrgbfolder = 'cleanspherergb';
cleandepthfolder = 'cleanspheredepth';
background_image_index = 1;
%good_image_indices = [1, 65, 138, 290, 333, 334, 335, 336 ...
good_image_indices = [65, 138, 336 ...
    417, 420, 421 ...
    ];

% good_image_indices = [69];

%Kernel ... ones indicate membership, zeros nonmembership.
% Centered about middle row and column
kernel_2D = ones(15,15);
num_samps = 2500;
scaling = .2;

assert(sum(sum(kernel_2D==0 | kernel_2D==1)) == numel(kernel_2D)); 

background_2D = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
background_2D = imresize(background_2D,scaling);
background_2D = double(background_2D)/255;

blur_halfwidth = size(background_2D,1)/1;
blur_filter = normpdf(floor(real(log(kron(exp(abs([-blur_halfwidth:blur_halfwidth])),exp(abs([-blur_halfwidth:blur_halfwidth]'))))) + .5),0,blur_halfwidth/4);
blur_filter = blur_filter.^5;

full_fig = [];

ref_imgs_2D = {};
recon_imgs_2D = {};

counter = 0;
for image_index=good_image_indices
    counter = counter + 1;
    
    depth_img_2D = imread([cleandepthfolder, '/ref',int2str(image_index),'.png']);
    depth_img_2D = double(depth_img_2D)/255;
    for color=1:3
        depth_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color),blur_filter,'same');
    end
    
    % Convert heightmap to r- and c- gradient maps.
    gradr_img_2D = conv2(depth_img_2D(:,:,1),[1;-1],'same');
    gradc_img_2D = conv2(depth_img_2D(:,:,1),[1 -1],'same');
    
    % Convert gradient maps to 3D normals maps
    normals_img_2D = zeros(size(depth_img_2D));
    normals_img_2D(:,:,1) = -gradr_img_2D;
    normals_img_2D(:,:,2) = -gradc_img_2D;
    normals_img_2D(:,:,3) = ones(size(depth_img_2D(:,:,3)));
    mags_img_2D = zeros(size(normals_img_2D(:,:,1)));
    for direction=1:3
        mags_img_2D = mags_img_2D + normals_img_2D(:,:,direction) .* ...
            normals_img_2D(:,:,direction);
    end
    for direction=1:3
        normals_img_2D(:,:,direction) = normals_img_2D(:,:,direction) ./ ...
            mags_img_2D;
    end
    
% % %     filter_2D_r_recon = imresize(filter_2D(:,:,:,1), size(depth_img_2D,1)/size(background_2D,1));
% % %     filter_2D_r_recon = filter_2D_r_recon / (size(depth_img_2D,1)/size(background_2D,1));
% % %     filter_2D_c_recon = imresize(filter_2D(:,:,:,2), size(depth_img_2D,1)/size(background_2D,1));
% % %     filter_2D_c_recon = filter_2D_c_recon / (size(depth_img_2D,1)/size(background_2D,1));
% % %     filter_2D_d_recon = imresize(filter_2D(:,:,:,3), size(depth_img_2D,1)/size(background_2D,1));
% % %     filter_2D_d_recon = filter_2D_d_recon / (size(depth_img_2D,1)/size(background_2D,1));
    
    ref_img_2D = imread([cleanrgbfolder, '/ref',int2str(image_index),'.png']);
    ref_img_2D = double(ref_img_2D)/255;

    %Reconstruct ref_img_2D from depth_img_2D
    recon_2D = zeros(size(depth_img_2D));
% % %     for color=1:3
% % %         recon_2D(:,:,color) = conv2(gradr_img_2D(:,:,color), filter_2D_r_recon(:,:,color), 'same') + ...
% % %             conv2(gradc_img_2D(:,:,color), filter_2D_c_recon(:,:,color), 'same') + ...
% % %             conv2(depth_img_2D(:,:,color), filter_2D_d_recon(:,:,color), 'same');
% % %     end
    for i=1:size(normals_img_2D,1)
        for j=1:size(normals_img_2D,2)
            for color=1:3
                norm_vec = zeros(3,1);
                norm_vec(:,1) = normals_img_2D(i,j,:);

                recon_2D(i,j,color) = full_model_lin(color,:) * norm_vec;
            end
        end
    end

    recon_2D = max(recon_2D, 0*recon_2D);
    background_2D_recon = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
    background_2D_recon = imresize(background_2D_recon, size(depth_img_2D,1)/size(background_2D_recon,1));
    background_2D_recon = double(background_2D_recon)/255;
%     recon_2D = recon_2D + background_2D_recon;

    grads_img_2D = zeros(size(depth_img_2D));
    grads_img_2D(:,:,1) = max(gradr_img_2D(:,:,1), 0*gradr_img_2D(:,:,1));
    grads_img_2D(:,:,2) = max(gradc_img_2D(:,:,1), 0*gradc_img_2D(:,:,1));
    
    next_fig = cat(1, grads_img_2D, ref_img_2D - background_2D_recon, recon_2D);
    if isempty(full_fig)
        full_fig = next_fig;
    else
        full_fig = cat(2, full_fig, next_fig);
    end
    
    ref_imgs_2D{counter} = ref_img_2D;
    recon_imgs_2D{counter} = recon_2D;
end


B = 0*eye(length(good_image_indices));
for i=1:length(good_image_indices)
    for j=1:length(good_image_indices)
        ref_img_2D = ref_imgs_2D{i};
        ref_img_2D = ref_img_2D - background_2D_recon;
        recon_img_2D = recon_imgs_2D{j};
        recon_img_2D = recon_img_2D - background_2D_recon;

        ref_img = reshape(ref_img_2D, [size(ref_img_2D,1)*size(ref_img_2D,2), 3]);
        recon_img = reshape(recon_img_2D, [size(recon_img_2D,1)*size(recon_img_2D,2), 3]);
        
        err_vec = zeros(3,1);
        for color=1:3
            err_vec(color,1) = ref_img(:,color)' * recon_img(:,color) / norm(ref_img) / norm(recon_img);
        end
        B(i,j) = norm(err_vec);
    end
end

disp(B);

figure;
imshow(full_fig);

filter = full_model_lin;
bg = background_2D;

end