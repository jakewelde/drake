function [filter, bg] = runGSFitBGSubNew()
% runGSFitBGSub uses the provided depthmaps and corresponding rgb images
% to train a linear convolution filter to map from one to the other.
%
% The linear model is as follows:
%
%   r/g/b_pixel = r/g/b_filter_2D' * depth_pixels
%
% Where depth_pixels is a patch of pixels from the source image and
% there are separate filters for r, g, and b.
%
% @retval filter the flipped 2D filter image, with separate color channels
% for r,g, and b stored in its third dimension. Can be applied to a depth
% image to produce one color=1/2/3 of an rgb image by means of the MATLAB
% command conv2, eg:
%
%    gs_image(:,:,color) = conv2(depth_image, filter(:,:,color), 'same');
%
% @retval bg the background that was subtracted from the images. Helpful
% for visualizing the derived GelSight image from a depth map, also shows
% what size of image the filter is meant to be applied to.

% cleanrgbfolder = 'cleanrgb';
% cleandepthfolder = 'cleandepth';
% background_image_index = 1;
% good_image_indices = [1, ...
%     5, ...
%     21, ...
%     69 ...
%     ];

figure;

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
num_samps = 2500+1250;
scaling = .2;

assert(sum(sum(kernel_2D==0 | kernel_2D==1)) == numel(kernel_2D)); 

out_imgs_2D = cell([1,length(good_image_indices)]);
in_imgs_2D = cell([1,length(good_image_indices)]);

background_2D = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
background_2D = imresize(background_2D,scaling);
background_2D = double(background_2D)/255;

blur_halfwidth = size(background_2D,1)/10;
gauss_blur_filter = normpdf(floor(real(log(kron(exp(abs([-blur_halfwidth:blur_halfwidth])),exp(abs([-blur_halfwidth:blur_halfwidth]'))))) + .5),0,blur_halfwidth/4);
gs_blur_filter = gauss_blur_filter .^ 4;

for image_counter=1:length(good_image_indices)
    true_image_index = good_image_indices(image_counter);
    
    ref_img_2D = imread([cleanrgbfolder, '/ref',int2str(true_image_index),'.png']);
    ref_img_2D = imresize(ref_img_2D,scaling);
    ref_img_2D = double(ref_img_2D)/255;
    ref_img_2D = ref_img_2D - background_2D;
    for color=1:3
        ref_img_2D(:,:,color) = conv2(ref_img_2D(:,:,color),fspecial('gaussian',floor(blur_halfwidth/4),blur_halfwidth/16),'same');
    end
    depth_img_2D = imread([cleandepthfolder, '/ref',int2str(true_image_index),'.png']);
    depth_img_2D = imresize(depth_img_2D,scaling);
    depth_img_2D = double(depth_img_2D)/255;
    for color=1:3
        depth_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color),gs_blur_filter,'same');
    end
    
    % Convert heightmap to r- and c- gradient maps.
    gradr_img_2D = conv2(depth_img_2D(:,:,1),[1;0;-1],'same');
    gradc_img_2D = conv2(depth_img_2D(:,:,1),[1,0,-1],'same');
    
    % Convert gradient maps to 3D normals maps
    normals_img_2D = zeros(size(depth_img_2D));
    normals_img_2D(:,:,1) = -gradr_img_2D;
    normals_img_2D(:,:,2) = -gradc_img_2D;
    normals_img_2D(:,:,3) = ones(size(depth_img_2D(:,:,3)));
    mags_img_2D = zeros(size(normals_img_2D(:,:,1)));
    for direction=1:3
        mags_img_2D = mags_img_2D + (normals_img_2D(:,:,direction) .* ...
            normals_img_2D(:,:,direction));
    end
    mags_img_2D = sqrt(mags_img_2D);
    for direction=1:3
        normals_img_2D(:,:,direction) = normals_img_2D(:,:,direction) ./ ...
            mags_img_2D;
    end
    
%     gradr_img_2D = conv2(depth_img_2D(:,:,1),[1,-1]','same');
%     gradc_img_2D = conv2(depth_img_2D(:,:,1),[1,-1],'same');
    
    disp(size(ref_img_2D));
    disp(length(out_imgs_2D));
    
    in_imgs_2D{image_counter} = ref_img_2D;
    out_imgs_2D{image_counter} = normals_img_2D;
%     
%     if image_counter==5
%         filter=depth_img_2D
%         return;
%     end
    
    imshow(ref_img_2D);
end

a = load('iamtrained.mat');
filter_2D = a.filter;

% filter_2D = trainConvLS(kernel_2D, in_imgs_2D, out_imgs_2D, num_samps);

full_fig = [];

for image_counter=1:length(good_image_indices)

    in_img_2D = in_imgs_2D{image_counter};
    out_img_2D = out_imgs_2D{image_counter};
 
    [heightmap_recon_2D, normals_recon_2D] = convertGStoHM(in_img_2D,filter_2D);
    
    
    next_fig = cat(1, in_img_2D, prettify_normal_img(out_img_2D), ...
        prettify_normal_img(-normals_recon_2D), repmat(heightmap_recon_2D,[1 1 3]));
    if isempty(full_fig)
        full_fig = next_fig;
    else
        full_fig = cat(2, full_fig, next_fig);
    end
    
    out_imgs_2D{image_counter} = out_img_2D;
end

figure;
imshow(full_fig);

filter = filter_2D;
bg = background_2D;

    function easy_to_see_2D=prettify_normal_img(normals_img_2D)
        easy_to_see_2D = normals_img_2D;
        easy_to_see_2D(:,:,3) = easy_to_see_2D(:,:,3) * 0;
        easy_to_see_2D = easy_to_see_2D / max(max(max(max(easy_to_see_2D))));
        %easy_to_see_2D(:,:,3) = abs(easy_to_see_2D(:,:,3) - 1);
        %easy_to_see_2D = min(1, easy_to_see_2D / max(max(max(max(easy_to_see_2D)))));
    end

end