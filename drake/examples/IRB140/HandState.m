classdef HandState < SingletonCoordinateFrame
  
  methods
    function obj=HandState(r, ind, name)
      typecheck(r,{'TimeSteppingRigidBodyManipulator','RigidBodyManipulator'});
      % Inds indicate which element of the overall state frame
      % our hand state is. 
      manipStateFrame = r.getManipulator().getStateFrame();
      manipStateFrame = manipStateFrame.getFrameByNum(ind);
      coordinates = manipStateFrame.getCoordinateNames();
      obj = obj@SingletonCoordinateFrame(name,length(coordinates),'x',coordinates);

      positionFrame = r.getManipulator().getPositionFrame();
      if getNumFrames(positionFrame)==1 && isempty(findTransform(obj,positionFrame))
        obj.addProjectionTransformByCoordinateNames(positionFrame);
      end

    end
  end
end
