function [c, dc] = staticStableEstimatorBox_SDFObjective(q, p, scanPts)
    % for every point in the scan, calculate SDF. Error is 
    % total quadratic error of SDF (they should all be zero in a perfect
    % fit)
    
    kinsol = p.doKinematics(q);
    [phi,n,x,x_body,body_idx] = p.signedDistances(kinsol,scanPts,false);
    
    K = 10 / size(scanPts, 2);
    c = K*phi.'*phi;
   % dc = 2*repmat(phi,1,3)*n;
    % derivatives
    bodies = unique(body_idx);
    dc = zeros(1,p.getNumPositions);
    for body=bodies.'
        [~, J] = p.forwardKin(kinsol, body, x_body);
        for i = find(body_idx == body).'
            dc = dc - K * 2 * phi(i) * n(:,i).'*J((i-1)*3+1:i*3, :);
        end
    end
%    lcmgl = LCMGLClient('sse_SDFOBJ_debug');
%    mindist = -0.1;
%    maxdist = 0.1;
%    phi(phi < mindist) = mindist;
%    phi(phi > maxdist) = maxdist;
%    norm_phi = (phi - mindist) / (maxdist - mindist);
%    for i = 1:size(phi)
%        lcmgl.glColor3f(1-norm_phi(i)^2,norm_phi(i)^2,0);
%        lcmgl.points(scanPts(1,i), scanPts(2,i), scanPts(3,i));
%    end
%    lcmgl.switchBuffers();

end