function [filter, bg] = runGSFitBGSubGradient()
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
num_samps = 5500;
scaling = .2;

assert(sum(sum(kernel_2D==0 | kernel_2D==1)) == numel(kernel_2D)); 

As = {};
bs = {};
ref_imgs = cell(length(good_image_indices),1);
depth_imgs = cell(length(good_image_indices),1);

background_2D = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
background_2D = imresize(background_2D,scaling);
background_2D = double(background_2D)/255;
    
for image_index=good_image_indices
    % Allocate the regression matrix A
    A = zeros(num_samps * length(good_image_indices),numel(kernel_2D),3);
    b = zeros(num_samps * length(good_image_indices),3);
    
    % Load in training images; heightmap + GelSight image
    ref_img_2D = imread([cleanrgbfolder, '/ref',int2str(image_index),'.png']);
    ref_img_2D = imresize(ref_img_2D,scaling);
    ref_img_2D = double(ref_img_2D)/255;
    depth_img_2D = imread([cleandepthfolder, '/ref',int2str(image_index),'.png']);
    depth_img_2D = imresize(depth_img_2D,scaling);
    depth_img_2D = double(depth_img_2D)/255;
    
    % Pre-processing; perform bg sub on GelSight image and
    % convert heightmap to r- and c- gradient maps.
    ref_img_2D = ref_img_2D - background_2D;
    gradr_img_2D = zeros(size(depth_img_2D));
    for color=1:3
        gradr_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color),[1,-1]','same');
    end
    gradc_img_2D = zeros(size(depth_img_2D));
    for color=1:3
        gradc_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color),[1,-1],'same');
    end
    
    disp(size(ref_img_2D))
    disp(size(depth_img_2D))
    
    assert(sum(size(ref_img_2D)~=size(depth_img_2D))==0);
    
    rows = size(ref_img_2D,1);
    cols = size(ref_img_2D,2);
    
    %Make 1D version of images and store them for later
    ref_img = reshape(ref_img_2D, [rows*cols, 3]);
    depth_img = reshape(depth_img_2D, [rows*cols, 3]);
    
    ref_imgs{image_index} = ref_img;
    depth_imgs{image_index} = depth_img;
    
    %Make linear kernel vector corresponding to starting kernel off in
    % top-left of image.
    krows = size(kernel_2D, 1);
    kcols = size(kernel_2D, 2);
    khfrows = floor(krows/2) + 1;
    khfcols = floor(kcols/2) + 1;
    kernel_lin = zeros(size(ref_img, 1),1); %zeros(krows*cols,1);
    for kcol=1:kcols
        kernel_lin(((kcol-1)*rows+ 1) : (kcol*rows)) = vertcat(kernel_2D(kcols,:)',zeros(rows-krows,1));
    end
    
    disp(size(kernel_lin)); % DEBUG
    disp(size(ref_img));
    disp(size(depth_img));
    
    %Randomly sample num_samps points from the image to build regression
    % matrix (will learn filter via least squares)
    sampd = zeros(num_samps,2);
    for i=1:num_samps
        roff = floor(rand()*(1+rows-krows));
        coff = floor(rand()*(1+cols-kcols));
        offset = coff*rows + roff;
        
        kernel_slid = vertcat(zeros(offset,1), kernel_lin(1:length(kernel_lin)-offset));
        
        r = khfrows + roff;
        c = khfcols + coff;
        
        in = depth_img(kernel_slid==1, :);
        out = ref_img_2D(r,c,:);
        
        A(i,:,:) = in;
        b(i,:) = out;
        
        sampd(i,:)=[r c];
    end
    disp(numel(in));
    disp(numel(out));
    %scatter(sampd(:,1),sampd(:,2));

    As{image_index}=A;
    bs{image_index}=b;
    
    imshow(ref_img_2D);
    
end

%Construct full matrix for regression
Aconv = vertcat(As{:});
bconv = vertcat(bs{:});

Afull = [];
bfull = [];
for color=1:3
    Afull(:,:,color) = blkdiag(Aconv(:,:,color)); %nothing else in this model
    bfull(:,color) = vertcat(bconv(:,color));  %nothing else in this model
end

disp(size(Afull));
disp(size(bfull));
disp(size(Aconv));
disp(size(bconv));

full_model_lin = zeros(size(Afull,2),3);
for color=1:3
    Amat = Afull(:,:,color);
    bvec = bfull(:,color);

    disp(size(Amat));
    disp(size(bvec));
    
    full_model_lin(:,color) = Amat \ bvec; %lsqlin(Amat, bvec, Amat*0, bvec*0);
end

filter_lin = full_model_lin(1:(krows*kcols),:,:);
filter_2D = reshape(filter_lin, krows, kcols, 3);
filter_2D = imrotate(filter_2D, 180);

%test_imgs = {'cleandepth/ref5.png','cleandepth/ref21.png','cleandepth/ref70.png'};

full_fig = [];

%for i=1:length(test_imgs)
%    img_name = test_imgs{i};
ref_imgs_2D = {};
recon_imgs_2D = {};
counter = 0;
for image_index=good_image_indices
    counter = counter + 1;
    depth_img_2D = imread([cleandepthfolder, '/ref',int2str(image_index),'.png']);
    filter_2D_recon = imresize(filter_2D, size(depth_img_2D,1)/size(background_2D,1));
    filter_2D_recon = filter_2D_recon / (size(depth_img_2D,1)/size(background_2D,1));% / (size(depth_img_2D,1)/size(background_2D,1));
    depth_img_2D = double(depth_img_2D)/255;
%     for color=1:3
%         depth_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color), [0,1,0;1,0,-1;0,-1,0], 'same');
%     end
    
    ref_img_2D = imread([cleanrgbfolder, '/ref',int2str(image_index),'.png']);
    ref_img_2D = double(ref_img_2D)/255;

    %Reconstruct ref_img_2D from depth_img_2D
    recon_2D = [];
    for color=1:3
        recon_color = conv2(depth_img_2D(:,:,color), filter_2D_recon(:,:,color), 'same');

        recon_2D = cat(3, recon_2D, recon_color);
    end

    %recon_2D = (recon_2D - min(min(min(recon_2D)))) / (max(max(max(recon_2D))) - min(min(min(recon_2D))));
    recon_2D = max(recon_2D, 0*recon_2D);
    background_2D_recon = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
    background_2D_recon = imresize(background_2D_recon, size(depth_img_2D,1)/size(background_2D_recon,1));
    background_2D_recon = double(background_2D_recon)/255;
    recon_2D = recon_2D + background_2D_recon;

    %subplot(1,3,1);
    %imshow(depth_img_2D);
    %subplot(1,3,2);
    %imshow(ref_img_2D);
    %subplot(1,3,3);
    %imshow(recon_2D);
    next_fig = cat(1, depth_img_2D, ref_img_2D, recon_2D);
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

figure;
imshow(filter_2D);%imresize(filter_2D,20));

filter = filter_2D;
bg = background_2D;

end