%function [filter, bg] = runGSFitBGSubNormals()
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

figure;

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

As = {};
bs = {};
ref_imgs = cell(length(good_image_indices),1);
depth_imgs = cell(length(good_image_indices),1);

background_2D = imread([cleanrgbfolder, '/ref',int2str(background_image_index),'.png']);
background_2D = imresize(background_2D,scaling);
background_2D = double(background_2D)/255;

blur_halfwidth = size(background_2D,1)/1;
blur_filter = normpdf(floor(real(log(kron(exp(abs([-blur_halfwidth:blur_halfwidth])),exp(abs([-blur_halfwidth:blur_halfwidth]'))))) + .5),0,blur_halfwidth/4);
blur_filter = blur_filter.^3;
blur_filter = blur_filter / sum(sum(sum(blur_filter)));

for image_index=good_image_indices
    % Allocate the regression matrix A
% % %     A = zeros(num_samps * length(good_image_indices),3*numel(kernel_2D),3);
% % %     b = zeros(num_samps * length(good_image_indices),3); %IS THIS RIGHT??
    A = zeros(size(background_2D,1)*size(background_2D,2),3,3);
    b = zeros(size(background_2D,1)*size(background_2D,2),3);
    
    % Load in training images; heightmap + GelSight image
    ref_img_2D = imread([cleanrgbfolder, '/ref',int2str(image_index),'.png']);
    ref_img_2D = imresize(ref_img_2D,scaling);
    ref_img_2D = double(ref_img_2D)/255;
    depth_img_2D = imread([cleandepthfolder, '/ref',int2str(image_index),'.png']);
    depth_img_2D = double(depth_img_2D)/255;
    for color=1:3
        depth_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color),blur_filter,'same');
    end
    depth_img_2D = imresize(depth_img_2D,scaling);
    
    % Pre-processing; perform bg sub on GelSight image
    ref_img_2D = ref_img_2D - background_2D;
    
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
        mags_img_2D = mags_img_2D + (normals_img_2D(:,:,direction) .* ...
            normals_img_2D(:,:,direction));
    end
    mags_img_2D = sqrt(mags_img_2D);
    for direction=1:3
        normals_img_2D(:,:,direction) = normals_img_2D(:,:,direction) ./ ...
            mags_img_2D;
    end
    
    assert(sum(sum(sum(normals_img_2D < -1))) == 0);
    assert(sum(sum(sum(normals_img_2D > 1))) == 0);
    for i=1:size(normals_img_2D,1)
        for j=1:size(normals_img_2D,2)
            vec = ones(3,1);
            vec(:,1) = normals_img_2D(i,j,:);
            mag=norm(vec);
            assert(.9 < mag && mag < 1.1, '%f and %f', mag, mags_img_2D(i,j));
        end
    end
    
    imshow((normals_img_2D + 1)/2);
    drawnow;
    
    disp(size(ref_img_2D))
    disp(size(depth_img_2D))
    
    assert(sum(size(ref_img_2D)~=size(depth_img_2D))==0);
    
    rows = size(ref_img_2D,1);
    cols = size(ref_img_2D,2);
    
    %Make 1D version of images and store them for later
    ref_img = reshape(ref_img_2D, [rows*cols, 3]);
    depth_img = reshape(depth_img_2D, [rows*cols, 3]);
%     gradr_img = reshape(gradr_img_2D, [rows*cols, 3]);
%     gradc_img = reshape(gradc_img_2D, [rows*cols, 3]);
    
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
% % %     sampd = zeros(num_samps,2);
% % %     for i=1:num_samps
% % %         roff = floor(rand()*(1+rows-krows));
% % %         coff = floor(rand()*(1+cols-kcols));
% % %         offset = coff*rows + roff;
% % %         
% % %         kernel_slid = vertcat(zeros(offset,1), kernel_lin(1:length(kernel_lin)-offset));
% % %         
% % %         r = khfrows + roff;
% % %         c = khfcols + coff;
% % %         
% % %         in = cat(1,gradr_img(kernel_slid==1,:),gradc_img(kernel_slid==1,:),depth_img(kernel_slid==1,:));
% % %         out = ref_img_2D(r,c,:);
% % %         
% % %         A(i,:,:) = in;
% % %         b(i,:) = out;
% % %         
% % %         sampd(i,:)=[r c];
% % %     end
    for i=1:size(normals_img_2D,1)
        for j=1:size(normals_img_2D,2)
            if sum(abs(background_2D(i,j,:)),3)>140/255
                norm_vec = zeros(3,1);
                norm_vec(:,1) = normals_img_2D(i,j,:);

                color_vec = zeros(3,1);
                color_vec(:,1) = ref_img_2D(i,j,:);

                A((i*size(background_2D,2))+j,:,:) = repmat(norm_vec,[1 3]);
                b((i*size(background_2D,2))+j,:) = color_vec;
            end
        end
    end

%     disp(numel(in));
%     disp(numel(out));

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

% % % filter_r_lin = full_model_lin(1:numel(kernel_2D),:,:);
% % % filter_r_2D = reshape(filter_r_lin, krows, kcols, 3);
% % % filter_r_2D = imrotate(filter_r_2D, 180);
% % % filter_c_lin = full_model_lin((numel(kernel_2D)+1):(2*numel(kernel_2D)),:,:);
% % % filter_c_2D = reshape(filter_c_lin, krows, kcols, 3);
% % % filter_c_2D = imrotate(filter_c_2D, 180);
% % % filter_d_lin = full_model_lin((2*numel(kernel_2D)+1):(3*numel(kernel_2D)),:,:);
% % % filter_d_2D = reshape(filter_d_lin, krows, kcols, 3);
% % % filter_d_2D = imrotate(filter_d_2D, 180);

% % % filter_2D = cat(4,filter_r_2D,filter_c_2D,filter_d_2D);

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

                recon_2D(i,j,:) = full_model_lin' * norm_vec;
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
    
    
    % Reconstruct the heightmap from the color image
    ref_img_bgsub = ref_img_2D - background_2D_recon;
    ref_img_bgsub = imresize(ref_img_bgsub,scaling);
    normals_recon_img = zeros(size(ref_img_bgsub));
    for i=1:size(normals_recon_img,1)
        for j=1:size(normals_recon_img,2)
            colorvec = zeros(3,1);
            colorvec(:,1) = ref_img_bgsub(i,j,:);
            normals_recon_img(i,j,:) = (full_model_lin')^(-1) * colorvec;
        end
    end
    mags_img_2D = zeros(size(normals_recon_img(:,:,1)));
    for direction=1:3
        mags_img_2D = mags_img_2D + (normals_recon_img(:,:,direction) .* ...
            normals_recon_img(:,:,direction));
    end
    mags_img_2D = sqrt(mags_img_2D);
    for direction=1:3
        normals_recon_img(:,:,direction) = normals_recon_img(:,:,direction) ./ ...
            mags_img_2D;
    end
    
    
    next_fig = cat(1, grads_img_2D, ref_img_2D - background_2D_recon, recon_2D);
    if isempty(full_fig)
        full_fig = next_fig;
    else
        full_fig = cat(2, full_fig, next_fig);
    end
    
    ref_imgs_2D{counter} = ref_img_2D;
    recon_imgs_2D{counter} = recon_img_2D;
    normals_imgs_2D{counter} = normals_recon_img;
end

% 
% B = 0*eye(length(good_image_indices));
% for i=1:length(good_image_indices)
%     for j=1:length(good_image_indices)
%         ref_img_2D = ref_imgs_2D{i};
%         ref_img_2D = ref_img_2D - background_2D_recon;
%         recon_img_2D = recon_imgs_2D{j};
%         recon_img_2D = recon_img_2D - background_2D_recon;
% 
%         ref_img = reshape(ref_img_2D, [size(ref_img_2D,1)*size(ref_img_2D,2), 3]);
%         recon_img = reshape(recon_img_2D, [size(recon_img_2D,1)*size(recon_img_2D,2), 3]);
%         
%         err_vec = zeros(3,1);
%         for color=1:3
%             err_vec(color,1) = ref_img(:,color)' * recon_img(:,color) / norm(ref_img) / norm(recon_img);
%         end
%         B(i,j) = norm(err_vec);
%     end
% end
% 
% disp(B);

figure;
imshow(full_fig);

figure;
one_img = normals_imgs_2D{1};
quiver(one_img(:,:,2),one_img(:,:,1));

filter = full_model_lin;
bg = background_2D;

%end