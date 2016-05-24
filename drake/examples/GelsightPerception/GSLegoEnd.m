classdef GSLegoEnd < RigidBodyManipulator
  
  methods
    
    function obj = GSLegoEnd()
        options = struct();
        options.floating = true;
        options.inertial = true;
        options.visual = true;
        obj = obj@RigidBodyManipulator(getFullPathFromRelativePath('urdf/legobrickend.urdf'), options);
        obj = addFrame(obj,RigidBodyFrame(findLinkId(obj,'base_link'),[0;0;0.0001],[0;0;0],'box_frame'));
        obj = compile(obj);
    end
    
  end
  
end