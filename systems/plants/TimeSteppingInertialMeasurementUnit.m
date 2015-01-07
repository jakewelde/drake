classdef TimeSteppingInertialMeasurementUnit < TimeSteppingRigidBodySensorWithState
  % IMU model that calculates acceleration numerically
  % from states over time
  
  methods
    function obj = TimeSteppingInertialMeasurementUnit(manip,body,xyz,rpy,noise)
      typecheck(body,'double');  % must be a body index
      
      if isa(manip,'PlanarRigidBodyManipulator')
        error('need to update kinematics implementation for planar manipulators (bug 1728) before this will work'); 
        if (nargin<3) xyz=zeros(2,1);
        else sizecheck(xyz,[2,1]); end
        if (nargin<4) rpy=0;
        else sizecheck(rpy,1); end
      else
        if (nargin<3) xyz=zeros(3,1);
        else sizecheck(xyz,[3,1]); end
        if (nargin<4) rpy=zeros(3,1);
        else sizecheck(rpy,[3,1]); end
        if (nargin<5) noise=0.0; end
      end      

      if (any(rpy)) error('Drake:TimeSteppingInertialMeasurementUnit:RPYNotImplemented', 'non-zero rpy not supported yet'); end  % but it wouldn't be too hard (see ContactForceTorque)
      
      obj.body = body;
      obj.xyz = xyz;
      obj.noise = noise;
      
      warning('the IMU outputs have not been tested yet (due to bug 1728)'); 
    end  

    function x0 = getInitialState(obj,tsmanip)
      warning('no real way to get init state?');
      x0 = zeros(24, 1);
    end
    
    function [obj,xdn,df] = update(obj,tsmanip,t,x,u)
      xdn = zeros(24, 1);
      xdn(13:24) = x(13:24);
      xdn(1:12) = x(1:12);
      df = [];
    end
    
    function y = output(obj,tsmanip,i,t,x,u)
      manip = tsmanip.getManipulator();
      this_x = x(13:24);
      last_x = x(25:end);
      
      % Calculate proper imu output
      numdof = getNumPositions(manip);
      q = this_x(1:numdof);
      qd = this_x(numdof+1:2*numdof);
      
      kinsol = doKinematics(manip,q,false,true,qd);
      [x,J] = forwardKin(manip,kinsol,obj.body,obj.xyz,1);
      Jdot = forwardJacDot(manip,kinsol,obj.body,obj.xyz,0,0);

      quat_body_to_world = rpy2quat(x(4:6));
      quat_world_to_body = quatConjugate(quat_body_to_world);
        
      % Get gyro output -- which we can get directly & accurately from
      % state
      rpy = x(4:6);
      rpydot = J(4:6,:)*qd;
      omega_base = rpydot2angularvel(rpy, rpydot); % TODO: replace with computation based on kinsol.twists
      omega_body = quatRotateVec(quat_world_to_body, omega_base);

      % Compute accelerometer in global frame from single timestep
      % numerically
      accel_base = (this_x(7:9) - last_x(7:9)) / tsmanip.timestep;
        accel_body = quatRotateVec(quat_world_to_body, accel_base);
        
      y = [omega_body; ...
        accel_body ];
      y = y + randn(size(y))*obj.noise/2;
    end
    
    function fr = constructFrame(obj,tsmanip)
      manip = tsmanip.getManipulator();
      body = getBody(manip,obj.body);
      fr = CoordinateFrame([strtok(body.linkname,'+'),'IMU'],6,'y', ...
        {'wx','wy','wz', ...       % angular rate (omega)
        'ax','ay','az'});         % linear acceleration
    end
    
    function fr = constructStateFrame(obj,tsmanip)
      manip = tsmanip.getManipulator();
      body = getBody(manip,obj.body);
      fr = CoordinateFrame([strtok(body.linkname,'+'),'OldState'],24,'y', ...
          {'x','y','z','r','p','y', ... 
          'xdot','ydot','zdot','rdot','pdot','ydot', ...
          'last_x','last_y','last_z','last_r','last_p','last_y', ... 
          'last_xdot','last_ydot','last_zdot','last_rdot','last_pdot','last_ydot'});
    end
    
    function tf = isDirectFeedthrough(obj)
      tf=false;
    end

    function obj = updateBodyIndices(obj,map_from_old_to_new)
      obj.body = map_from_old_to_new(obj.body);
    end
  end
  
  properties
    body
    xyz
    noise
  end
  
end