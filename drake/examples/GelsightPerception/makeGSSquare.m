function [gs_manifold_x, gs_manifold_n]=makeGSSquare(w, width, gs_thickness)
    % Make a standard rectangular GS sensor, xyz-centered on the origin
    if nargin < length({'w', 'width', 'gs_thickness'})
        gs_thickness = 0.05;
    end
    if nargin < length({'w', 'width'})
        width=1.0;
    end
    if nargin < length({'w'})
        w=200;
    end
    
    row = width * ( floor(-w/2)+1:floor(w/2) ) / w; % a row [-w, ..., w]
    grid = cat(1,repmat(row,1,w),kron(row,ones(1,w))); % a 2-by-w*w mat, row x row:
     % [[-w, ..., w,   -w, ..., w,  .......  , -w, ..., w]
     %  [-w, ...,-w, -w+1,...,-w+1, .......  ,  w, ..., w]]
    gs_manifold_x = cat(1,grid,zeros(1,w*w)); % 3-by-w*w, col-vecs of positions
    gs_manifold_n = repmat([0;0;gs_thickness], 1, w*w); % 3-by-w*w, col-vecs of norms
end