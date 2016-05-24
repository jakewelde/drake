classdef GSSimpleBox < RigidBodyManipulator
  
  methods
    
    function obj = GSSimpleBox()
        options = struct();
        options.floating = true;
        options.inertial = true;
        options.visual = true;
        obj = obj@RigidBodyManipulator(getFullPathFromRelativePath('urdf/justabox.urdf'), options);
        obj = addFrame(obj,RigidBodyFrame(findLinkId(obj,'base_link'),[0;0;.1],[0;0;0],'box_frame'));
        obj = compile(obj);
    end
    
  end
  
end