function [gs_manifold_x, gs_manifold_n]=makeGSRectangular(pixelw, pixelh, width, height, thickness)
    % Make a standard rectangular GS sensor, xyz-centered on the origin
    if nargin < length({'w', 'width', 'gs_thickness'})
        thickness = 0.05;
    end
    if nargin < length({'w', 'width'})
        width=1.0;
    end
    if nargin < length({'w'})
        pixelw=200;
    end
    
    row_h = height * ( floor(-pixelh/2)+1:floor(pixelh/2) ) / pixelh; % a row [-h, ..., h]*height/h
    row_w = width * ( floor(-pixelw/2)+1:floor(pixelw/2) ) / pixelw; % a row [-w, ..., w]*width/h
    grid = cat(1,repmat(row_h,1,pixelw),kron(row_w,ones(1,pixelh))); % a 2-by-w*h mat, row x row:
     % [[-w, ..., w,   -w, ..., w,  .......  , -w, ..., w] .* [[width/w ]
     %  [-h, ...,-h, -h+1,...,-h+1, .......  ,  h, ..., h]] .* [height/w]]
    gs_manifold_x = cat(1,grid,zeros(1,pixelh*pixelw)); % 3-by-w*h, col-vecs of positions
    gs_manifold_n = repmat([0;0;thickness], 1, pixelh*pixelw); % 3-by-w*h, col-vecs of norms
end