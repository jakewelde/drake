function prob=probMatchGelSight(depth_img_2D, measured_img_2D, gs_bg, gs_filter)
% collisionGelSight  Determines the probability that the given depth
%   map corresponds with the given measured GelSight image. Background
%   and filter are required for model to work.
%
% @param depth_img_2D Size height by width by 3.
% @param measured_img_2D Size height by width by 3.
% @param gs_bg Size height by width by 3.
% @param gs_filter Size n by m by 3. Should be such that convolving
% depth_img_2D by gs_filter and adding gs_bg gives a decent approximation
% to the true GelSight image for gs_filter.
%
% @retval prob probability that the two images correspond.

% Produce backgroundless predicted image from depth image
pred_img_2D = zeros(size(depth_img_2D,1), size(depth_img_2D,2), 1);
for color=1:3
    pred_img_2D(:,:,color) = conv2(depth_img_2D(:,:,color), gs_filter(:,:,color), 'same');
end

% Compute score as the cosine-distance between the predicted and measured
% image
measured_img_2D = measured_img_2D - gs_bg;
pred_img = reshape(pred_img_2D, [size(pred_img_2D,1)*size(pred_img_2D,2), 3]);
measured_img = reshape(measured_img_2D, [size(measured_img_2D,1)*size(measured_img_2D,2), 3]);

corr_vec = zeros(3,1);
for color=1:3
    corr_vec(color,1) = pred_img(:,color)' * measured_img(:,color) / norm(pred_img) / norm(measured_img);
end
score = norm(corr_vec); %max(0,.01+sum(corr_vec));

% Convert this score into a probability and return it
prob = score + .0001;

if isnan(prob)
    prob = .0001;
end

end