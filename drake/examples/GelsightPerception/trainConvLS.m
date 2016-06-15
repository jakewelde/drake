function [filter] = trainConvLS(kernel_mask, input_imgs, output_imgs, samps_per_img)
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

Arows = samps_per_img * size(output_imgs{1},3);
Acols = numel(kernel_mask) * size(input_imgs{1},3) * size(output_imgs{1},3);

Aconv = zeros(Arows * length(input_imgs),Acols);
bconv = zeros(Arows * length(input_imgs),1);

for image_index=1:length(input_imgs)
    in_img_2D = input_imgs{image_index};
    out_img_2D = output_imgs{image_index};

    % allocate the regression matrix A, which will do r,g,b-patch ->
    % r,g,b-pixel
    A = zeros(Arows, Acols);
    b = zeros(Arows,1);
    
    assert(sum(size(out_img_2D(:,:,1))~=size(in_img_2D(:,:,1)))==0, 'Images at index %d not same-sized!', image_index);
    
    rows = size(out_img_2D,1);
    cols = size(out_img_2D,2);
    
    %Make 1D version of images
    out_img = reshape(out_img_2D, [rows*cols, size(out_img_2D,3)]);
    in_img = reshape(in_img_2D, [rows*cols, size(in_img_2D,3)]);
    
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
        
        in = zeros(numel(kernel_mask),size(in_img_2D,3));
        out = zeros(size(out_img_2D,3),1);
        in(:,:) = in_img(kernel_slid==1, :);
        out(:) = out_img_2D(r,c,:);
        
        % in (size k-by-3?) now stores a numel(kernel_2D)-sized patch of r,g,b values from the input image,
        % and
        % out (size 1) now stores a single r,g,b value from the output image. We can
        % use these to build our regression.
        Achunk = zeros(size(out_img_2D,3), Acols + size(in_img_2D,3)*size(out_img_2D,3)); %some padding for easier assign
        Arow = zeros(1, Acols);
        for j=1:size(in_img_2D,3)
            ej = zeros(1,size(in_img_2D,3)*size(out_img_2D,3));
            ej(j) = 1;
            Arow = Arow + kron(in(:,j)',ej);
        end
        for j=1:size(out_img_2D,3)
            horiz_off = (j-1)*size(in_img_2D,3)+1;
            Achunk(j,(horiz_off):(horiz_off + Acols - 1)) = Arow;
        end
        Achunk = Achunk(:,1:Acols);
        
        rownum = ((i-1)*size(out_img_2D,3) + 1);
        
        A((rownum) : (rownum+size(out_img_2D,3)-1),:) = Achunk;
        b((rownum) : (rownum+size(out_img_2D,3)-1)) = out(:);
        
        sampd(i,:)=[r c];
    end
    disp(numel(in));
    disp(numel(out));
%     
%     disp(A(1:12,1:12));
%     disp(b(1:12));

    Aconv((Arows*(image_index-1) + 1):(Arows*(image_index-1) + Arows),:) = A;
    bconv((Arows*(image_index-1) + 1):(Arows*(image_index-1) + Arows)) = b;
    
    imshow(out_img_2D);
    
end

%Keep some images around, for size reference
in_img_2D = input_imgs{image_index};
out_img_2D = output_imgs{image_index};

% Afull = zeros(size(Aconv));
% bfull = zeros(size(bconv));
% 
% Afull(:,:) = blkdiag(Aconv(:,:)); %nothing else in this model
% bfull(:) = vertcat(bconv(:));  %nothing else in this model
% 
% disp(size(Afull));
% disp(size(bfull));
% disp(size(Aconv));
% disp(size(bconv));
% 
% disp(size(Afull));
% disp(size(bfull));

full_model_vec = Aconv \ bconv;
% full_model_vec is the whole, ugly vector. Want to reshape into
% convolution kernels indexed by in_color,out_color.

full_model_lin = reshape(full_model_vec, [size(in_img_2D,3), size(out_img_2D,3), numel(kernel_mask)]);
full_model_lin = shiftdim(full_model_lin,2);
% Can index as full_model_lin(pixel_in_kernel, in_color, out_color);

full_model_2D = zeros(size(kernel_mask,1),size(kernel_mask,2),size(in_img_2D,3),size(out_img_2D,3));
for in_color = 1:size(in_img_2D,3)
    for out_color = 1:size(out_img_2D,3)
        one_model_lin = full_model_lin(:,in_color,out_color);
        full_model_2D(:,:,in_color,out_color) = reshape(one_model_lin,[size(kernel_mask,1),size(kernel_mask,2)]);
        full_model_2D(:,:,in_color,out_color) = imrotate(full_model_2D(:,:,in_color,out_color),180);
    end
end
% Finally, can index a 2D conv kernel as full_model_2D(:,:,in_color,out_color)

disp(size(full_model_2D));
disp(size(full_model_lin));
disp(size(kernel_mask));
disp(size(in_img_2D,3));
disp(size(out_img_2D,3));

filter = full_model_2D;

end