classdef BallSLAMThing < DrakeSystem
  
  methods 
    function obj = BallSLAMThing(ballplant)
      typecheck(ballplant,'BallPlant2D');
      obj = obj@DrakeSystem(0, 0, 2, 2);  % number of inputs = 2
      
      sys = BallFlightPhasePlant2D();
      obj = setInputFrame(obj,sys.getOutputFrame);
      obj = setOutputFrame(obj,sys.getOutputFrame);
      
    end
    
  end
  
  properties
  end
end
