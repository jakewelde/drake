function runGSFit

refs_list = dir('reference');

% good_image_indices = [1,2,3,4, ...
%     5,6,7,8, ...
%     21,22,23,24, ...
%     69,70,71,72, ...
%     ];

good_image_indices = [1, ...
    5, ...
    21, ...
    69 ...
    ];

% good_image_indices = [69];

%Kernel ... ones indicate membership, zeros nonmembership.
% Centered about middle row and column
kernel_2D = ones(10,10);
num_samps = 12000;

assert(sum(sum(kernel_2D==0 | kernel_2D==1)) == numel(kernel_2D)); 

% allocate the regression matrix A
A = zeros(num_samps * length(good_image_indices),numel(kernel_2D),3);
b = zeros(num_samps * length(good_image_indices),1,3);

for image_index=good_image_indices
    ref_img_2D = imread(['cleanrgb/ref',int2str(image_index),'.png']);
    depth_img_2D = imread(['cleandepth/ref',int2str(image_index),'.png']);
    ref_img_2D = imresize(ref_img_2D,.2);
    depth_img_2D = imresize(depth_img_2D,.2);
    
    ref_img_2D = double(ref_img_2D)/255;
    depth_img_2D = double(depth_img_2D)/255;
    
    disp(size(ref_img_2D))
    disp(size(depth_img_2D))
    
    assert(sum(size(ref_img_2D)~=size(depth_img_2D))==0);
    
    rows = size(ref_img_2D,1);
    cols = size(ref_img_2D,2);
    
    %Make 1D version of images
    ref_img = reshape(ref_img_2D, [rows*cols, 3]);
    depth_img = reshape(depth_img_2D, [rows*cols, 3]);
    
    %Make kernel vector corresponding to starting kernel off in top-left
    % of image.
    krows = size(kernel_2D, 1);
    kcols = size(kernel_2D, 2);
    khfrows = floor(krows/2) + 1;
    khfcols = floor(kcols/2) + 1;
    
%     rows = 6;
%     cols = 3;
%     ref_img = zeros(rows*cols, 3);
%     depth_img = zeros(rows*cols, 3);
    
    kernel_lin = zeros(krows*cols,1);
    for kcol=1:kcols
        kernel_lin(((kcol-1)*rows+ 1) : (kcol*rows)) = vertcat(kernel_2D(kcols,:)',zeros(rows-krows,1));
    end
    kernel_lin = vertcat(kernel_lin, zeros(size(ref_img, 1) - size(kernel_lin, 1),1));
    
    disp(size(kernel_lin));
    
    disp(size(ref_img));
    disp(size(depth_img));
    
    %Randomly sample num_samps points from the image
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
        b(i,1,:) = out;
        
        %figure;
        %imshow(reshape(in, krows, kcols, 3));
                
        sampd(i,:)=[r c];
    end
    disp(numel(in));
    disp(numel(out));
    %scatter(sampd(:,1),sampd(:,2));
    
    
    filter_lin = zeros(krows*kcols,3);
    for color=1:3
        Amat = A(:,:,color);
        bvec = b(:,:,color);
        
        filter_lin(:,color) = lsqlin(Amat, bvec, Amat*0, bvec*0); %Amat \ bvec;
    end
        
    filter_2D = reshape(filter_lin, krows, kcols, 3);
    
    %imshow(filter_2D); 
    
    %Reconstruct ref_img_2D from depth_img_2D
    recon_2D = [];
    for color=1:3
        recon_color = conv2(depth_img_2D(:,:,color), filter_2D(:,:,color));
        
        recon_2D = cat(3, recon_2D, recon_color);
    end
    recon_2D = recon_2D(krows:krows+rows-1, kcols:kcols+cols-1, :);
    
    %recon_2D = (recon_2D - min(min(min(recon_2D)))) / (max(max(max(recon_2D))) - min(min(min(recon_2D))));
    recon_2D = max(recon_2D, 0*recon_2D);
    
    figure
    %subplot(1,3,1);
    %imshow(depth_img_2D);
    %subplot(1,3,2);
    %imshow(ref_img_2D);
    %subplot(1,3,3);
    %imshow(recon_2D);
    imshow(cat(1, depth_img_2D, ref_img_2D, recon_2D));
    
% 	figure;
%   plot(filter_lin);
    
    %ref_img_2D = (ref_img_2D(:,:,1) + ref_img_2D(:,:,2) + ref_img_2D(:,:,3))/3;
    %depth_img_2D = (depth_img_2D(:,:,1) + depth_img_2D(:,:,2) + depth_img_2D(:,:,3))/3;
    
    %ref_img = padarray(ref_img, [10 10], 0, 'both');
    
    %imshow(filter);
end

end