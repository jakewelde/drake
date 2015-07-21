classdef RigidBodySlopedTerrain < RigidBodyTerrain
%  This class provides an implementation of RigidBodyTerrain with z increasing
%  in the x direction, with z specified at the origin, and slope specified
  
  methods 
    function obj = RigidBodySlopedTerrain(varargin)
      % Construct a sloped terrain map
      p = inputParser();
      p.addOptional('z', 0, @isscalar);
      p.addOptional('dzdx', 0.1, @isscalar);
      p.parse(varargin{:});
      obj.z = p.Results.z;
      obj.dzdx = p.Results.dzdx;
      obj.geom = constructRigidBodyGeometry(obj);
    end
    
    function [z,normal] = getHeight(obj,xy)
      [m n] = size(xy);
      z = repmat(obj.z,1,n);
      z = z + xy(1, :)*obj.dzdx;
      normal = repmat([-obj.dzdx;0;1],1,n);
    end

    function geom = getCollisionGeometry(obj)
      geom = obj.geom;
    end

    function geom = getVisualGeometry(obj)
      geom = obj.geom;
    end

    function obj = setGeometryColor(obj, color)
      geom = obj.getVisualGeometry();
      geom.c = reshape(color, 3, 1);
      obj.geom = geom;
    end
    
    function geom = constructRigidBodyGeometry(obj)
      box_width = 100;
      box_depth = 1;
      geom = RigidBodyBox([box_width;box_width;box_depth], [box_depth/2*sin(obj.dzdx),0,obj.z-(box_depth/2)*cos(obj.dzdx)], [0, -atan2(obj.dzdx, 1), 0]);
      geom.c = hex2dec({'ee','cb','ad'})'/256;  % something a little brighter (peach puff 2 from http://www.tayloredmktg.com/rgb/);
      geom.name = 'terrain';
%      geom.c = hex2dec({'cd','af','95'})'/256;
    end
  end
  
  properties
    geom;
    z;
    dzdx;
  end
end
