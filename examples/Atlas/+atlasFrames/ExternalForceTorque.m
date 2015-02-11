classdef ExternalForceTorque < LCMCoordinateFrame & Singleton
  
  methods
    function obj=ExternalForceTorque()

      coordinates = {'body_or_frame_id','fx','fy','fz','tx','ty','tz'};

      obj = obj@LCMCoordinateFrame('ExternalForceTorque',length(coordinates),'f');
      obj = obj@Singleton();
      
%       if isempty(obj.lcmcoder)
%         coder = drc.control.ForceTorqueStateCoder();
%         setLCMCoder(obj,JLCMCoder(coder));
%         obj.setCoordinateNames(coordinates);
%         obj.setDefaultChannel('EXTERNAL_FORCE_TORQUE');
%       end
    end
  end
end