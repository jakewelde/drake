classdef ExternalForceTorque < LCMCoordinateFrame & Singleton
  
  methods
    function obj=ExternalForceTorque()
      coordinates = {'body_or_frame_id','fx','fy','fz','tx','ty','tz'};
      obj = obj@LCMCoordinateFrame('ExternalForceTorque',drake.lcmt_external_force_torque,'f',coordinates);
      obj = obj@Singleton();
    end
  end
end