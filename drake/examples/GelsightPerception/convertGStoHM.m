function [heightmap_recon_2D, normals_recon_2D]=convertGStoHM(gs_img_2D, gs_filter, use_qp)
    % convertGStoHM transforms the provided rgb GelSight image into a heightmap
    % by making use of filter_2D, a k-by-k-by-3-by-3 convolution filter trained
    % to take a patch of rgb pixels to an xyz normal vector.
    %
    % @param gs_img_2D a background-subtracted RGB GelSight image, size [w h 3], with color values
    % normalized to the range [0,1].
    %
    % @param filter_2D RGB->NxNyNz convolution filter, size [k1 k2 3 3], which
    % takes a k1-by-k2-sized patch of cells from gs_img_2D and maps them to the
    % normals from which the heightmap is reconstructed.
    %
    % @param use_qp tells whether to use QP or least squares for solve.
    %
    % @retval heightmap_recon_2D the heightmap inferred from the image.
    %
    % @retval normals_recon_2D the unnormalized normal vectors produced by the
    % convolution of arg1 and arg2.
    
    if nargin < 3
        use_qp = true;
    end

    recon_2D = zeros(size(gs_img_2D));
    for color1=1:size(gs_img_2D,3)
        for color2=1:size(gs_img_2D,3)
            recon_2D(:,:,color2) = recon_2D(:,:,color2) + conv2(gs_img_2D(:,:,color1), gs_filter(:,:,color1,color2), 'same');
        end
    end

    hsize = floor(size(recon_2D,1)/10);
    sigma = floor(hsize/8);
    for color=1:3
        recon_2D(:,:,color) = conv2(recon_2D(:,:,color),fspecial('gaussian',hsize,sigma),'same');
    end

    %% Recover un-magnitude-corrected normals from convolution
    normals_recon_2D = recon_2D;

    gradr_recon_2D = zeros(size(normals_recon_2D,1),size(normals_recon_2D,2));
    gradr_recon_2D(:,:) = -normals_recon_2D(:,:,1);
    gradc_recon_2D = zeros(size(normals_recon_2D,1),size(normals_recon_2D,2));
    gradc_recon_2D(:,:) = -normals_recon_2D(:,:,2);

    grads_recon_2D = 100*cat(3,gradr_recon_2D,gradc_recon_2D,zeros(size(gradr_recon_2D))); % + .0094;

    %% Build a QP to find a nice heightmap

    gradr_vec = cat(1,zeros(1,size(gradr_recon_2D,2)),gradr_recon_2D,zeros(1,size(gradr_recon_2D,2)));
    gradr_vec = gradr_vec(:);
    gradc_vec = cat(2,zeros(size(gradc_recon_2D,1),1),gradc_recon_2D,zeros(size(gradc_recon_2D,1),1));
    gradc_vec = gradc_vec(:);

    %These sparse matrices encode the row- and col-wise gradient operations
    Dr = convmtx2([1;0;-1],[size(grads_recon_2D,1),size(grads_recon_2D,2)]);
    Dc = convmtx2([1,0,-1],[size(grads_recon_2D,1),size(grads_recon_2D,2)]);

    border_left_2D = cat(2,ones(size(grads_recon_2D,1),1),zeros(size(grads_recon_2D,1),size(grads_recon_2D,2)-1))==1;
    border_top_2D = cat(1,ones(1,size(grads_recon_2D,2)),zeros(size(grads_recon_2D,1)-1,size(grads_recon_2D,2)))==1;
    border_full_2D = (0*grads_recon_2D(:,:,1))==1; %initially full of false
    border_full_2D = border_full_2D | border_left_2D; %set one border
    border_full_2D = border_full_2D | imrotate(border_left_2D,180); %set one border
    border_full_2D = border_full_2D | border_top_2D; %set one border
    border_full_2D = border_full_2D | imrotate(border_top_2D,180); %set one border
    border_full_2D = double(border_full_2D);
    Q = spdiags(border_full_2D(:),0,numel(border_full_2D),numel(border_full_2D));

    if (use_qp)
        H = Dr'*Dr + Dc'*Dc + Q'*Q;
        f = -2*gradr_vec'*Dr + -2*gradc_vec'*Dc;

        x0 = 0*gradr_recon_2D(:);
        x0 = x0 + .001;

        opts = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','none');
        heightmap_recon = quadprog(H,f,[],[],[],[],[],[],[],opts);
    else
        heightmap_recon = cat(1,Dr,Dc,Q) \ cat(1,gradr_vec,gradc_vec,zeros(size(Q,2),1));
    end
    
    heightmap_recon_2D = reshape(heightmap_recon,size(gradr_recon_2D));
    heightmap_recon_2D = max(0,heightmap_recon_2D);
    heightmap_recon_2D = 300*heightmap_recon_2D;

end