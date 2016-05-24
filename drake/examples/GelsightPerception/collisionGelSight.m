function gs_depths=collisionGelSight(rbm, kinsol, gs_manifold_x, gs_manifold_n, use_margins)
% collisionGelSight  Produce simulated GelSight depth readings.
%   gs_depths=collisionGelSight(rbm, kinsol, gs_manifold_x, gs_manifold_n,
%   use_margins) produces a vector of depth readings, one for each point on
%   the gelsight manifold gs_manifold_x (with thickness vectors
%   gs_manifold_n) by casting rays onto the given RigidBodyManipulator.
%
%   See also collisionApproachGelSight,
%   RigidBodyManipulator/collisionRaycast.
%
% @param rbm a RigidBodyManipulator to collide against
%
% @param kinsol result of calling doKinematics(obj, q) where q is a
% position vector.
%
% @param gs_manifold_x position vectors corresponding to the unyielding
% "back" manifold of the GelSight sensor. Size is 3 x N.
%
% @param gs_manifold_n normal vectors corresponding to the direction and
% thickness of the GelSight sensor at every point in its manifold. Size is
% 3 x N.
%
% @param use_margins boolean indicating whether or not to use a collision
%   model whose boxes and meshes are padded with a margin to improve
%   numerical stability of contact gradient. Default true.
%
%
% @retval gs_depth distances reported by the simulated GelSight sensor.
% Size is N x 1.

if nargin < length({'rbm', 'kinsol', 'gs_manifold_x', 'gs_manifold_n','use_margins'})
   use_margins = true; 
end

gs_endpoints_x = gs_manifold_x + gs_manifold_n;
raw_depths = rbm.collisionRaycast(kinsol, gs_manifold_x, gs_endpoints_x, use_margins);

gs_depths = raw_depths;

gs_thicknesses = (gs_manifold_n(1,:).*gs_manifold_n(1,:));
gs_thicknesses = gs_thicknesses + (gs_manifold_n(2,:).*gs_manifold_n(2,:));
gs_thicknesses = gs_thicknesses + (gs_manifold_n(3,:).*gs_manifold_n(3,:));
gs_thicknesses = sqrt(gs_thicknesses');

gs_depths(gs_depths==-1) = gs_thicknesses(gs_depths==-1);

assert(sum(gs_depths > (gs_thicknesses + .01)) == 0);
assert(min(gs_depths) >= 0 - .01);

end