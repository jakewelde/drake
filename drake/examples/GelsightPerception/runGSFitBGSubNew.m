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

%filter_2D = trainConvLS(kernel_2D, in_imgs_2D, out_imgs_2D, num_samps);

full_fig = [];

for image_counter=1:length(good_image_indices)

    in_img_2D = in_imgs_2D{image_counter};
    out_img_2D = out_imgs_2D{image_counter};
    recon_2D = zeros(size(out_img_2D));
    for color1=1:size(in_img_2D,3)
        for color2=1:size(out_img_2D,3)
            recon_2D(:,:,color2) = recon_2D(:,:,color2) + conv2(in_img_2D(:,:,color1), filter_2D(:,:,color1,color2), 'same');
        end
    end
    
    hsize = floor(size(recon_2D,1)/10);
    sigma = floor(hsize/8);
    for color=1:3
        recon_2D(:,:,color) = conv2(recon_2D(:,:,color),fspecial('gaussian',hsize,sigma),'same');
    end
    %recon_2D = max(0, min(1, recon_2D));
    
    % Use normals to reconstruct a heightmap
    normals_recon_2D = recon_2D;
% % %     mags_img_2D = zeros(size(normals_recon_2D(:,:,1)));
% % %     for direction=1:3
% % %         mags_img_2D = mags_img_2D + (normals_recon_2D(:,:,direction) .* ...
% % %             normals_recon_2D(:,:,direction));
% % %     end
% % %     mags_img_2D = sqrt(mags_img_2D);
% % %     for direction=1:3
% % %         normals_recon_2D(:,:,direction) = normals_recon_2D(:,:,direction) ./ ...
% % %             mags_img_2D;
% % %     end
    
    assert(min(min(min(normals_recon_2D(:,:,3))))<=0, 'normal z comp 0 or less: %f',min(min(min(normals_recon_2D(:,:,3)))));
    assert(min(min(min(normals_recon_2D(:,:,3))))<=0.1, 'normal z comp 0.1 or less: %f',min(min(min(normals_recon_2D(:,:,3)))));
    
    gradr_recon_2D = zeros(size(normals_recon_2D,1),size(normals_recon_2D,2));
    gradr_recon_2D(:,:) = -normals_recon_2D(:,:,1);% ./ normals_recon_2D(:,:,3);
    gradc_recon_2D = zeros(size(normals_recon_2D,1),size(normals_recon_2D,2));
    gradc_recon_2D(:,:) = -normals_recon_2D(:,:,2);% ./ normals_recon_2D(:,:,3);
    
    
    grads_recon_2D = 100*cat(3,gradr_recon_2D,gradc_recon_2D,zeros(size(gradr_recon_2D))); % + .0094;
    
    %grads_recon_2D = max(0,min(0.1,grads_recon_2D));
    %grads_recon_2D = (grads_recon_2D-min(min(min(grads_recon_2D))))/(max(max(max(grads_recon_2D)))-min(min(min(grads_recon_2D))));
    
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
    'Algorithm','interior-point-convex','Display','iter');
    heightmap_recon = quadprog(H,f,[],[],[],[],[],[],[],opts);
    
    heightmap_recon_2D = reshape(heightmap_recon,size(gradr_img_2D));
    heightmap_recon_2D = heightmap_recon_2D*300;
%     heightmap_recon_2D = (heightmap_recon_2D-min(min(min(heightmap_recon_2D)))) / ...
%         (max(max(max(heightmap_recon_2D)))-min(min(min(heightmap_recon_2D))));
    disp('Dr:');
    disp(size(gradr_vec));
    disp(size(Dr));
    
    disp('Dc:');
    disp(size(gradc_vec));
    disp(size(Dc));
    
    disp(max(max(max(grads_recon_2D))));
    disp(min(min(min(grads_recon_2D))));
    disp(mean(mean(mean(grads_recon_2D))));
    disp(sum(sum(sum(grads_recon_2D > .0094)))/numel(grads_recon_2D));
    %next_fig = cat(1, repmat(in_img_2D,[1 1 3]), out_img_2D, grads_recon_2D);
    %next_fig = cat(1, in_img_2D, out_img_2D, grads_recon_2D);
    
    next_fig = cat(1, in_img_2D, prettify_normal_img(out_img_2D), ...
        grads_recon_2D, repmat(heightmap_recon_2D,[1 1 3]));
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