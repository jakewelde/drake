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
    
    ref_img_2D = ref_img_2D - background_2D;
    
    disp(size(ref_img_2D));
    disp(length(out_imgs_2D));
    
    in_imgs_2D{image_counter} = ref_img_2D;
    out_imgs_2D{image_counter} = normals_img_2D;
    
    imshow(ref_img_2D);
end

% a = load('iamtrained.mat');

%filter_2D = a.filter;
filter_2D = trainConvLS(kernel_2D, in_imgs_2D, out_imgs_2D, num_samps);

full_fig = [];

for image_counter=1:length(good_image_indices)
%     counter = counter + 1;
%     depth_img_2D = imread([cleandepthfolder, '/ref',int2str(true_image_index),'.png']);
%     filter_2D_recon = imresize(filter_2D, size(depth_img_2D,1)/size(background_2D,1));
%     filter_2D_recon = filter_2D_recon / (size(depth_img_2D,1)/size(background_2D,1));% / (size(depth_img_2D,1)/size(background_2D,1));
%     depth_img_2D = double(depth_img_2D)/255;
% %     for color=1:3
% %         depth_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color), [0,1,0;1,0,-1;0,-1,0], 'same');
% %     end
%     
%     ref_img_2D = imread([cleanrgbfolder, '/ref',int2str(true_image_index),'.png']);
%     ref_img_2D = double(ref_img_2D)/255;
% 
%     %Reconstruct ref_img_2D from depth_img_2D
%     recon_2D = [];
%     for color=1:3
%         recon_color = conv2(depth_img_2D(:,:,color), filter_2D_recon(:,:,color), 'same');
% 
%         recon_2D = cat(3, recon_2D, recon_color);
%     end
%     % recon_2D = recon_2D(krows:krows+rows-1, kcols:kcols+cols-1, :);
% 
%     %recon_2D = (recon_2D - min(min(min(recon_2D)))) / (max(max(max(recon_2D))) - min(min(min(recon_2D))));
%     recon_2D = max(recon_2D, 0*recon_2D);
%     background_2D_recon = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
%     background_2D_recon = imresize(background_2D_recon, size(depth_img_2D,1)/size(background_2D_recon,1));
%     background_2D_recon = double(background_2D_recon)/255;
%     recon_2D = recon_2D + background_2D_recon;

    in_img_2D = in_imgs_2D{image_counter};
    out_img_2D = out_imgs_2D{image_counter};
    recon_2D = zeros(size(out_img_2D));
    for color1=1:size(in_img_2D,3)
        for color2=1:size(out_img_2D,3)
            recon_2D(:,:,color2) = recon_2D(:,:,color2) + conv2(in_img_2D(:,:,color1), filter_2D(:,:,color1,color2), 'same');
        end
    end
    
    for color=1:3
        recon_2D(:,:,color) = conv2(recon_2D(:,:,color),gauss_blur_filter,'same');
    end
    recon_2D = max(0, min(1, recon_2D));
    
    %next_fig = cat(1, repmat(in_img_2D,[1 1 3]), out_img_2D, recon_2D);
    %next_fig = cat(1, in_img_2D, out_img_2D, recon_2D);
    
    next_fig = cat(1, in_img_2D, prettify_normal_img(out_img_2D), prettify_normal_img(recon_2D));
    if isempty(full_fig)
        full_fig = next_fig;
    else
        full_fig = cat(2, full_fig, next_fig);
    end
    
    out_imgs_2D{image_counter} = out_img_2D;
    recon_imgs_2D{image_counter} = recon_2D;
end

figure;
imshow(full_fig);

B = 0*eye(length(good_image_indices));
for i=1:length(good_image_indices)
    for j=1:length(good_image_indices)
        out_img_2D = out_imgs_2D{i};
        recon_img_2D = recon_imgs_2D{j};

        ref_img = reshape(out_img_2D, [size(out_img_2D,1)*size(out_img_2D,2), 3]);
        recon_img = reshape(recon_img_2D, [size(recon_img_2D,1)*size(recon_img_2D,2), 3]);
        
        err_vec = zeros(3,1);
        for color=1:3
            err_vec(color,1) = ref_img(:,color)' * recon_img(:,color) / norm(ref_img) / norm(recon_img);
        end
        B(i,j) = norm(err_vec);
    end
end

disp(B);

% figure;
% imshow(filter_2D);%imresize(filter_2D,20));

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