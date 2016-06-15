function prob=probMatchGelSightNaive(depth_img_2D, measured_img_2D, gs_bg)
% collisionGelSight  Determines the probability that the given depth
%   map corresponds with the given measured GelSight image. Background
%   and filter are required for model to work.
%
% @param depth_img_2D Size height by width by 3.
% @param measured_img_2D Size height by width by 3.
% @param gs_bg Size height by width by 3.
%
% @retval prob probability that the two images correspond.

% Produce backgroundless predicted image from depth image
pred_img_2D = zeros(size(depth_img_2D,1), size(depth_img_2D,2), 1);
gradr_img_2D = conv2(depth_img_2D(:,:,1),[1,-1]','same');
gradc_img_2D = conv2(depth_img_2D(:,:,1),[1,-1],'same');

normals_img_2D = zeros(size(depth_img_2D));
normals_img_2D(:,:,1) = -gradr_img_2D;
normals_img_2D(:,:,2) = -gradc_img_2D;
normals_img_2D(:,:,3) = ones(size(depth_img_2D(:,:,3)));
norms_img_2D = zeros(size(normals_img_2D(:,:,1)));
for direction=1:3
    norms_img_2D = norms_img_2D + normals_img_2D(:,:,direction) .* ...
        normals_img_2D(:,:,direction);
end
for direction=1:3
    normals_img_2D(:,:,direction) = normals_img_2D(:,:,direction) ./ ...
        norms_img_2D;
end

%Dot-product norms with 45-degree vectors
mags_img_2D = sin(2*acos(normals_img_2D(:,:,3)));

% Compute score as the cosine-distance between the predicted and measured
% image
measured_img_2D = measured_img_2D - gs_bg;
measured_img_2D = max(measured_img_2D, 0);

measured_int_img_2D = 0*measured_img_2D(:,:,1);
for color=1:3
    measured_int_img_2D = measured_img_2D(:,:,color).*measured_img_2D(:,:,color);
end
measured_int_img_2D = measured_int_img_2D - .00001;


% Blur for better gradients
blur_filter = fspecial('gaussian',floor([size(gs_bg,1)/10, size(gs_bg,2)/10]),size(gs_bg,1)/30);
mags_img_2D = conv2(mags_img_2D, blur_filter, 'same');
measured_int_img_2D = conv2(measured_int_img_2D, blur_filter, 'same');

% figure;
% imshow(mags_img_2D);
% figure;
% imshow(measured_int_img_2D);

mags_img = reshape(mags_img_2D, [size(mags_img_2D,1)*size(mags_img_2D,2),1]);
measured_int_img = reshape(measured_int_img_2D, [size(measured_int_img_2D,1)*size(measured_int_img_2D,2),1]);

score = mags_img' * measured_int_img / norm(mags_img) / norm(measured_int_img);

% Convert this score into a probability and return it
%prob = score + .0001;
prob = 2^score;

if isnan(prob)
    prob = .0001;
end

end