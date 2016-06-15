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
num_samps = 250;
scaling = .2;

assert(sum(sum(kernel_2D==0 | kernel_2D==1)) == numel(kernel_2D)); 

ref_imgs_2D = cell([1,length(good_image_indices)]);
depth_imgs_2D = cell([1,length(good_image_indices)]);

background_2D = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
background_2D = imresize(background_2D,scaling);
background_2D = double(background_2D)/255;

for image_counter=1:length(good_image_indices)
    true_image_index = good_image_indices(image_counter);
    
    ref_img_2D = imread([cleanrgbfolder, '/ref',int2str(true_image_index),'.png']);
    ref_img_2D = imresize(ref_img_2D,scaling);
    ref_img_2D = double(ref_img_2D)/255;
    depth_img_2D = imread([cleandepthfolder, '/ref',int2str(true_image_index),'.png']);
    depth_img_2D = imresize(depth_img_2D,scaling);
    depth_img_2D = double(depth_img_2D)/255;
    
    gradr_img_2D = conv2(depth_img_2D(:,:,1),[1,-1]','same');
    gradc_img_2D = conv2(depth_img_2D(:,:,1),[1,-1],'same');
    for color=1:3
        depth_img_2D(:,:,color) = (gradr_img_2D .* gradr_img_2D) + (gradc_img_2D .* gradc_img_2D);
    end
    
    ref_img_2D = ref_img_2D - background_2D;
    
    disp(size(ref_img_2D));
    disp(length(ref_imgs_2D));
    
    ref_imgs_2D{image_counter} = ref_img_2D;
    depth_imgs_2D{image_counter} = depth_img_2D;
    
    imshow(ref_img_2D);
end

filter_2D = trainConvLS(kernel_2D, ref_imgs_2D, depth_imgs_2D, num_samps);


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

    depth_img_2D = depth_imgs_2D{image_counter};
    ref_img_2D = ref_imgs_2D{image_counter};
    recon_2D = zeros(size(ref_img_2D));
    for color=1:3
        for color2=1:3
            recon_2D(:,:,color) = recon_2D(:,:,color) + conv2(depth_img_2D(:,:,color2), filter_2D(:,:,color2), 'same');
        end
        recon_2D(:,:,color) = recon_2D(:,:,color) / 3;
    end
    
    next_fig = cat(1, depth_img_2D, ref_img_2D, recon_2D);
    if isempty(full_fig)
        full_fig = next_fig;
    else
        full_fig = cat(2, full_fig, next_fig);
    end
    
    ref_imgs_2D{image_counter} = ref_img_2D;
    recon_imgs_2D{image_counter} = recon_2D;
end


B = 0*eye(length(good_image_indices));
for i=1:length(good_image_indices)
    for j=1:length(good_image_indices)
        ref_img_2D = ref_imgs_2D{i};
        recon_img_2D = recon_imgs_2D{j};

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

figure;
imshow(filter_2D);%imresize(filter_2D,20));

filter = filter_2D;
bg = background_2D;

end