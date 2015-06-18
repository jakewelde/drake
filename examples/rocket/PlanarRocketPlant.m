classdef PlanarRocketPlant < DrakeSystem

  % state:  
  %  q(1) - x position
  %  q(2) - z position
  %  q(3) - pitch (theta)
  % input:
  %  u(1) - thrust
  %  u(2) - gimbal

  properties
    T = [0; -2.0]; % vector from com to thrust application pt
    ldraw = 4; % rocket length, for drawing purposes
    wdraw = 0.5; % rocket width, for drawing purposes
    m = 1.0; % mass of rocket
    I = 1.0;
    g = 9.81; % gravity
    Kt = 1.0; % thrust coefficient
    
    num_q = 3; % x, y, theta
    % todo: drag?
    
    thrust_noise_profile;
    gimbal_noise_profile;
    noise_profile_w = [0.5, 1.0, 2.0, 4.0, 8.0];
    
    x0 = [0.0; 15.0; 0.0; 0.0; 0.0; 0];
    
    noise = 0;
  end
  
  methods
    function obj = PlanarRocketPlant()
      obj = obj@DrakeSystem(6,0,2,6,false,true);
      obj = obj.setOutputFrame(obj.getStateFrame);  % allow full-state feedback
      
      obj.thrust_noise_profile = (rand([length(obj.noise_profile_w), 1]) - 0.5)*2;
      obj.gimbal_noise_profile = (rand([length(obj.noise_profile_w), 1]) - 0.5)*2;
      
    end
    
    function [xdot, dxdot] = dynamics(obj,t,x,u)

      if (obj.noise)
        u_noise = [0;0];
        for i=1:length(obj.noise_profile_w)
          u_noise(1) = u_noise(1) + cos(t*obj.noise_profile_w(i))*0.5;
          u_noise(2) = u_noise(2) + cos(t*obj.noise_profile_w(i))*0.05;
        end
        u = u + u_noise;
      end
      
      %[x(1:3) [u(1:2);0]]
      
      q=x(1:obj.num_q); 
      qd=x((obj.num_q+1):end);
      % thrust is applied along y axis plus gimbal amount
      thrust = u(1);
      gimbal = u(2);
      
      thrust_vector = [obj.Kt * thrust * sin(q(3)+gimbal);
                       obj.Kt * thrust * cos(q(3)+gimbal)];
      torque = cross([rotmat(-q(3))*obj.T;0], [thrust_vector;0]);
      torque = torque(3);
      qdd = [ thrust_vector / obj.m - [0; obj.g];
        torque / obj.I];
      
      xdot = [qd;qdd];
      
      if (nargout==2)
        % let's run out some derivatives...
        % no dependence on x/y, but dependence on theta, thrust, and gimbal
        dthrust_vector_dtheta = [obj.Kt * thrust * cos(q(3)+gimbal);
                                 -obj.Kt * thrust * sin(q(3)+gimbal)];
        dthetadd_dtheta = (-sin(q(3))*obj.T(1) + cos(q(3))*obj.T(2))*thrust_vector(2) + ...
                          (cos(q(3))*obj.T(1) + sin(q(3))*obj.T(2))*dthrust_vector_dtheta(2) + ...
                          (-cos(q(3))*obj.T(1) + sin(q(3))*obj.T(2))*thrust_vector(1) + ...
                          (-sin(q(3))*obj.T(1) - cos(q(3))*obj.T(2))*dthrust_vector_dtheta(1);
        dqdd_dtheta = [dthrust_vector_dtheta / obj.m; dthetadd_dtheta / obj.I];
        if (thrust == 0)
          dqdd_dthrust = zeros(3, 1);
        else
          dqdd_dthrust = [thrust_vector / thrust / obj.m; 
                          ((cos(q(3))*obj.T(1) + sin(q(3))*obj.T(2))*thrust_vector(2)/thrust - ...
                          (sin(q(3))*obj.T(1) + cos(q(3))*obj.T(2))*thrust_vector(1)/thrust)/ obj.I];
        end
        dqdd_dgimbal = [dthrust_vector_dtheta / obj.m;
                        ((cos(q(3))*obj.T(1) + sin(q(3))*obj.T(2))*dthrust_vector_dtheta(2) - ...
                        (sin(q(3))*obj.T(1) + cos(q(3))*obj.T(2))*dthrust_vector_dtheta(1))/ obj.I];
        
        dxdot = [zeros(obj.num_q,1+obj.num_q), eye(obj.num_q), zeros(obj.num_q,obj.num_u);
                 zeros(obj.num_q,1+2), dqdd_dtheta, zeros(obj.num_q, 3), dqdd_dthrust, dqdd_dgimbal];
               
        % numerical:
%         eps = 1E-6;
%         dxdot_num = zeros(obj.num_q*2, 1+obj.num_q*2+obj.num_u);
%         for i=1:obj.num_q*2
%           diff = zeros(obj.num_q*2, 1);
%           diff(i) = eps;
%           dxdot_num(:, i+1) = (dynamics(obj,t,x+diff,u) - xdot) / eps;
%         end
%         for i=1:obj.num_u
%           diff = zeros(obj.num_u, 1);
%           diff(i) = eps;
%           dxdot_num(:, obj.num_q*2+i+1) = (dynamics(obj,t,x,u+diff) - xdot) / eps;
%         end    
%         
%         err = norm(dxdot_num - dxdot)
      end
      
      
    end
    
    function y = output(obj,t,x,u)
      % default output is the full state
      y = x;
    end
    
    function x = getInitialState(obj)
      x = obj.x0;
    end
    
    function obj = setInitialState(obj, x0)
      obj.x0 = x0;
    end
    
    function [c,V] = hoverLQR(obj, xdes)
      xdes = Point(obj.getStateFrame,xdes);
      u0 = Point(obj.getInputFrame,obj.m*obj.g/obj.Kt * [1;0]);
      Q = diag([10 10 10 0.1 0.1 0.1]);
      R = [10.0 0.00; 0.00 1.0];  %R = diag([0.1 0.1]);

      options.angle_flag = [0;0;1;0;0;0];
      if (nargout>1)
        [c,V0] = tilqr(obj,xdes,u0,Q,R);
        sys = feedback(obj,c);

        pp = sys.taylorApprox(0,xdes,[],3);  % make polynomial approximation
        options.degL1=2;
        %options.method='bilinear';
        %options.degV=4;
        V=regionOfAttraction(pp,V0,options);
      else
        c = tilqr(obj,xdes,u0,Q,R, options);
      end
    end
    
    function [c] = hoverTVLQR(obj, xtraj, T)
      utraj = setOutputFrame(PPTrajectory(zoh([0, T+1], obj.m*obj.g/obj.Kt * [1 1;0 0])), obj.getInputFrame);
      Q = diag([10 10 10 0.1 0.1 0.1]);
      Qf = Q;
      R = [10.0 0.00; 0.00 1.0];  %R = diag([0.1 0.1]);

      options.angle_flag = [0;0;1;0;0;0];
      c = tvlqr(obj,xtraj,utraj,Q,R, Qf, options);
    end
    
  end
  
end
