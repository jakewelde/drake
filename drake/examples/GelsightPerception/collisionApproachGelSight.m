function [distance,collided]=collisionApproachGelSight(rbm, kinsol, gs_manifold_x, approach_vec, use_margins)
% collisionApproachGelSight  Tells how far the gs manifold can be moved
% along the given approach vector before the GelSight sensor has
% maximal-penetration contact with the world.
%   distance, collided=collisionApproachGelSight(rbm, kinsol, gs_manifold_x,
%   approach_vec) produces distance corresponding to how far along
%   approach_vec the gelsight manifold gs_manifold_x can be moved before it
%   collides with the given RigidBodyManipulator. collided encodes the
%   index of the point that was hit.
%
%   See also collisionGelSight, RigidBodyManipulator/collisionRaycast.
%
% @param rbm a RigidBodyManipulator to collide against
%
% @param kinsol result of calling doKinematics(rbm, q) where q is a
% position vector.
%
% @param gs_manifold_x position vectors corresponding to the unyielding
% "back" manifold of the GelSight sensor. Size is 3 x N.
%
% @param approach_vec vector of size 3 indicating the path that the
% GelSight will be moved along until it collides with rbm. This is NOT a
% normalized vector; the sensor will be moved some amount along the LENGTH
% of approach_vec. Size is 3 x 1.
%
% @param use_margins boolean indicating whether or not to use a collision
%   model whose boxes and meshes are padded with a margin to improve
%   numerical stability of contact gradient. Default true.
%
%
% @retval distance min distance the sensor can be translated before it
% collides with rbm, or -1 if it will never collide.
%
% @retval collided index of the first manifold point to collide with the
% hull, or -1 if it will never collide.

if nargin < length({'rbm', 'kinsol', 'gs_manifold_x', 'approach_vec','use_margins'})
   use_margins = true; 
end

gs_endpoints_x = gs_manifold_x + repmat(approach_vec, 1, size(gs_manifold_x,2));
raw_depths = rbm.collisionRaycast(kinsol, gs_manifold_x, gs_endpoints_x, use_margins);

if (sum(raw_depths==-1) == 0)
    distance = -1;
    collided = -1;
    return
else
    raw_depths(raw_depths==-1) = norm(approach_vec)+.000001;
    [distance, collided] = min(raw_depths);
    return
end

end