function [filter] = trainGreyConvLS(kernel_mask, input_imgs, output_imgs, samps_per_img)
% trainConvLS trains a linear convolution filter to map from one set of input
% images to a corresponding set of output images.
%
% @retval filter the flipped 2D filter image, with separate color channels
% for r,g, and b stored in its third dimension. Can be applied to a novel 
% input image to produce one color=1/2/3 of an rgb output image by means of the MATLAB
% command conv2, eg:
%
%   out_image(:,:,color) = conv2(in_image, filter(:,:,color), 'same');

if nargin < length({'kernel_mask','input_imgs','output_imgs','samps_per_img'})
    samps_per_img = 2500;
end

assert(length(input_imgs)==length(output_imgs), ...
    'Input image cell and output image cell must be of same length.');
assert(sum(sum(kernel_mask==1))==numel(kernel_mask), ...
    'Mask must have entries equal only to 1.')

%Kernel ... ones indicate membership, zeros nonmembership.
% Centered about middle row and column
kernel_mask = ones(15,15);

As = cell([1,length(input_imgs)]);
bs = cell([1,length(input_imgs)]);

for image_index=1:length(input_imgs)
    % allocate the regression matrix A
    A = zeros(samps_per_img * length(input_imgs),numel(kernel_mask),3);
    b = zeros(samps_per_img * length(input_imgs),3);
    
    in_img_2D = input_imgs{image_index};
    out_img_2D = output_imgs{image_index};
    
    assert(sum(size(out_img_2D)~=size(in_img_2D))==0, 'Images at index %d not same-sized!', image_index);
    
    rows = size(out_img_2D,1);
    cols = size(out_img_2D,2);
    
    %Make 1D version of images
    out_img = reshape(out_img_2D, [rows*cols, 3]);
    in_img = reshape(in_img_2D, [rows*cols, 3]);
    
    %Make kernel vector corresponding to starting kernel off in top-left
    % of image.
    krows = size(kernel_mask, 1);
    kcols = size(kernel_mask, 2);
    khfrows = floor(krows/2) + 1;
    khfcols = floor(kcols/2) + 1;
    
    disp(image_index);
    disp(size(out_img));
    
    kernel_lin = zeros(size(out_img, 1),1);
    for kcol=1:kcols
        kernel_lin(((kcol-1)*rows+ 1) : (kcol*rows)) = vertcat(kernel_mask(kcols,:)',zeros(rows-krows,1));
    end
    
    disp(size(kernel_lin));
    
    disp(size(out_img));
    disp(size(in_img));
    
    %Randomly sample num_samps points from the image
    sampd = zeros(samps_per_img,2);
    for i=1:samps_per_img
        roff = floor(rand()*(1+rows-krows));
        coff = floor(rand()*(1+cols-kcols));
        offset = coff*rows + roff;
        
        kernel_slid = vertcat(zeros(offset,1), kernel_lin(1:length(kernel_lin)-offset));
        
        r = khfrows + roff;
        c = khfcols + coff;
        
        in = in_img(kernel_slid==1, :);
        out = out_img_2D(r,c,:);
        
        A(i,:,:) = in;
        b(i,:) = out;
        
        sampd(i,:)=[r c];
    end
    disp(numel(in));
    disp(numel(out));

    As{image_index}=A;
    bs{image_index}=b;
    
    imshow(out_img_2D);
    
end

%Construct full matrix for regression
Aconv = vertcat(As{:});
bconv = vertcat(bs{:});

Afull = zeros(size(Aconv));
bfull = zeros(size(bconv));
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
    
    full_model_lin(:,color) = Amat \ bvec;
end

filter_lin = full_model_lin(1:(krows*kcols),:,:);
filter_2D = reshape(filter_lin, krows, kcols, 3);
filter_2D = imrotate(filter_2D, 180);

filter = filter_2D;

end