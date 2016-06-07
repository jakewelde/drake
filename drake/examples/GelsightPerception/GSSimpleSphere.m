classdef GSSimpleSphere < RigidBodyManipulator
  
  methods
    
    function obj = GSSimpleSphere()
        options = struct();
        options.floating = true;
        options.inertial = true;
        options.visual = true;
        obj = obj@RigidBodyManipulator(getFullPathFromRelativePath('urdf/justasphere.urdf'), options);
        obj = addFrame(obj,RigidBodyFrame(findLinkId(obj,'base_link'),[0;0;0],[0;0;0],'main_frame'));
        obj = compile(obj);
    end
    
  end
  
end