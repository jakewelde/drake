classdef Rocket < RigidBodyManipulator
  properties
      lcmgl = LCMGLClient('rocket');
  end
    
  methods
    
    function obj = Rocket(options)
      w = warning('off','Drake:RigidBodyManipulator:ReplacedCylinder');
      %warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
      obj = obj@RigidBodyManipulator(getFullPathFromRelativePath('rocket.urdf'),options);
      warning(w);
      
      obj = addSensor(obj,FullStateFeedbackSensor);
      obj = compile(obj);
    end
   
    function I = getInertia(obj)
      I = obj.body(2).inertia;
    end
    
    function u0 = nominalThrust(obj)
      % thrust command m*g
      u0 = Point(getInputFrame(obj),[0;0;getMass(obj)*norm(getGravity(obj))*ones(1,1); 0.1]);
    end
    
    function traj_opt = addPlanVisualizer(obj,traj_opt)
      % spew out an lcmgl visualization of the trajectory.  intended to be
      % used as a callback (fake objective) in the direct trajectory
      % optimization classes

      if ~checkDependency('lcmgl')
        warning('lcmgl dependency is missing.  skipping visualization'); 
        return;wrench_lcmgl = LCMGLClient('wrench');
      end
      lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'RocketPlan');
      
      typecheck(traj_opt,'DirectTrajectoryOptimization');

      traj_opt = traj_opt.addDisplayFunction(@(x)visualizePlan(x,lcmgl),traj_opt.x_inds(1:3,:));
      
      function visualizePlan(x,lcmgl)
        lcmgl.glColor3f(1, 0, 0);
        lcmgl.glPointSize(3);
        lcmgl.points(x(1,:),x(2,:),x(3,:));
        lcmgl.glColor3f(.5, .5, 1);
        lcmgl.plot3(x(1,:),x(2,:),x(3,:));
        lcmgl.switchBuffers;
      end
    end
    
    function drawThrust(obj, t, inp)
      u = inp(1:4);
      x = inp(5:end);

      kinsol = doKinematics(obj, x(1:8));
      start_pos = obj.forwardKin(kinsol, 4, [0;0;0], 0);
      arrow = obj.forwardKin(kinsol, 4, [0;0;-u(3)*0.1], 0) - start_pos;
      
      head_width = 0.05;
      head_length = 0.05;
      body_width = 0.03;
      color = [1.0, 0.0, 0.0];
      
      length = norm(arrow);
      obj.lcmgl.glLineWidth(2);
      obj.lcmgl.glPushMatrix();
      obj.lcmgl.glTranslated(start_pos(1),start_pos(2),start_pos(3));
      obj.lcmgl.glColor3f(color(1),color(2),color(3));
%             lcmgl.glBegin(lcmgl.LCMGL_LINES);
%             lcmgl.glVertex3f(0,0,0);
%             lcmgl.glVertex3f(force_ij(1)/obj.force_scaler,force_ij(2)/obj.force_scaler,force_ij(3)/obj.force_scaler);
%             lcmgl.glEnd();
      % compute the rotation angle
      rotate_axis = cross([1;0;0],arrow/length);
      if(norm(rotate_axis)>0.01)
        rotate_angle = asin(norm(rotate_axis))/pi*180;
        rotate_axis = rotate_axis/norm(rotate_axis);
        obj.lcmgl.glRotated(rotate_angle,rotate_axis(1),rotate_axis(2),rotate_axis(3));
      end
      
      obj.lcmgl.glPushMatrix();
      obj.lcmgl.glRotated(90,0,1,0);
      obj.lcmgl.cylinder([0;0;0],body_width,body_width,length-head_length,20,20);
      obj.lcmgl.glTranslated(0,0,length-head_length);
      obj.lcmgl.cylinder([0;0;0],head_width,0,head_length,20,20);
      obj.lcmgl.glPopMatrix();
%             lcmgl.drawArrow3d(norm(force_ij)/obj.force_scaler,0.03,0.03,0.01);
      obj.lcmgl.glPopMatrix();
      obj.lcmgl.switchBuffers();
    end
  end
  
  
  methods (Static)
    
    
    function runLQR
      options.floating = true;
      options.terrain = RigidBodyFlatTerrain();
      r = Rocket(options); 
      sys = TimeSteppingRigidBodyManipulator(r,.001, options);
      
      v = sys.constructVisualizer();
      v_thr = FunctionHandleVisualizer(MultiCoordinateFrame({r.getInputFrame, r.getOutputFrame.frame{1}, r.getOutputFrame.frame{2}}), @r.drawThrust);
      x0 = [0;0;5.0;zeros(13,1)];
      u0 = double(r.nominalThrust());
      [A,B] = linearize(r,0,x0,u0);
      
      Q = diag([1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0 ...
                1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0]);
      R = diag([1. 1. 0.1 0.1]);
      K = lqr(full(A),full(B),Q,R);
     
      c = AffineSystem([],[],[],[],[],[],[],-K,u0 + K*x0);
      c = setInputFrame(c,getStateFrame(r));
      c = setOutputFrame(c,getInputFrame(r));
 
      outs(1).system = 1;
      outs(1).output = 1;
      outs(2).system = 1;
      outs(2).output = 2;
      outs(3).system = 2;
      outs(3).output = 1;
      sys = mimoFeedback(sys,c, [], [], [], outs);
      
      outs(1).system = 1;
      outs(1).output = 1;
      outs(2).system = 1;
      outs(2).output = 2;
      outs(3).system = 1;
      outs(3).output = 3;
      
      sys = mimoCascade(sys, v_thr, [], [], outs);
      sys = mimoCascade(sys,v, [], [], outs);
      
      x0(4:6) = x0(4:6) + 0.1*rand(3, 1);
      x0(3) = x0(3) + 3;
      x0(1) = x0(1) + 5;
      simulate(sys,[0 6],double(x0));
    end
    function runOpenLoop
      options.floating = true;
      options.terrain = RigidBodyFlatTerrain();
      r = Rocket(options); 
      sys = TimeSteppingRigidBodyManipulator(r,.001, options);
      
      v = sys.constructVisualizer();

      x0 = [0;0;4.0;zeros(13,1)];
      u0 = Point(sys.getInputFrame, [0.0; 0.0; 5.0, 0.0]);
      
      sys = cascade(ConstantTrajectory(u0),sys);

      sys = cascade(sys,v);
      x0(4:6) = x0(4:6) + 0.1*rand(3, 1);
      simulate(sys,[0 10],double(x0));
      
      %options.capture_lcm_channels = 'LCMGL';
      %[ytraj,xtraj,lcmlog] = simulate(sys,[0 0.05],double(x0),options);
      %lcmlog
      %v.playback(xtraj,struct('lcmlog',lcmlog));
%      figure(1); clf; fnplt(ytraj);
    end
  end
end
