classdef GSSimpleBox < RigidBodyManipulator
  
  methods
    
    function obj = GSSimpleBox()
        options = struct();
        obj = obj@RigidBodyManipulator(getFullPathFromRelativePath('justabox.urdf'), options);
        obj = compile(obj);
    end
    
  end
  
end