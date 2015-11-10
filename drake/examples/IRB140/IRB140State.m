classdef IRB140State < SingletonCoordinateFrame
  
  methods
    function obj=IRB140State(r)
      typecheck(r,{'TimeSteppingRigidBodyManipulator','RigidBodyManipulator'});
      manipStateFrame = r.getManipulator().getStateFrame();
      if (isa(manipStateFrame.getFrameByNum(1), 'MultiCoordinateFrame'))
        manipStateFrame = manipStateFrame.getFrameByNum(1);
      end
      coordinates = manipStateFrame.getCoordinateNames();
      obj = obj@SingletonCoordinateFrame('IRB140State',length(coordinates),'x',coordinates);
      
      positionFrame = r.getManipulator().getPositionFrame();
      if getNumFrames(positionFrame)==1 && isempty(findTransform(obj,positionFrame))
        obj.addProjectionTransformByCoordinateNames(positionFrame);
      end
    end
  end
end
