classdef PFModeHopping < DrakeSystem
    
  properties
      x0; % initial state of the system we're estimating
      nq; % 
      nx; % state size of system we're estimating
      nu; % inputs of system we're estimating
      nz; % outputs of system we're estimating
      plant; % system we're estimating
      v; % viewer of it
      
      nC; % # of contact points
      nD; % # of friction cone members
      nContactForces;  
      
      lcmgl;
      
      num_particles;
      
      % STATE VECTOR:
      % [ particle_weights; % num_particles x 1, should be normalized
      %   particle_1_state; % (nX + nContactForces) x 1
      %   ...
      %   particle_<num_particles>_state; % (nX + nContactForces) x 1
      % ]
  end
  
  methods
    function obj = PFModeHopping(plant, num_particles, x0, state_frame)
      typecheck(plant,'TimeSteppingRigidBodyManipulator');
      
      nx = getNumStates(plant);
      nu = getNumInputs(plant);
      nz = getNumOutputs(plant);
      
      % contact info
      multiple_contacts = false;
      nC = plant.getNumContactPairs;
      [~,normal,d] = plant.contactConstraints(zeros(plant.getNumPositions,1), multiple_contacts);
      nD = 2*length(d);
      assert(size(normal,2) == nC); % just a double check
      % contact forces along friction cone bases
      nContactForces = nC*(1 + nD);
      
      %             weights             state, contact state
      num_states = num_particles + num_particles*(nx+nContactForces);
      obj = obj@DrakeSystem(0,num_states,nu+nz,num_states,0,1);
      obj.num_particles = num_particles;
      
      obj.x0 = x0;
      obj.plant = plant;
      obj.v = obj.plant.constructVisualizer();
      obj.nx = nx;
      obj.nq = obj.plant.getNumPositions();
      obj.nu = nu;
      obj.nz = nz;
      obj.nC = nC;
      obj.nD = nD;
      obj.nContactForces = nContactForces;
      
      % setup input/output frames
      plant_input_frames = getInputFrame(plant);
      plant_output_frames = getOutputFrame(plant);
      
      % this should be replaced with a generic flatten-multicoordframe
      if obj.nu > 0
          if isa(plant_input_frames, 'MultiCoordinateFrame')
            frames = plant_input_frames.frame;
          else
            frames = {plant_input_frames};
          end
      else
          frames = {};
      end
      if isa(plant_output_frames, 'MultiCoordinateFrame')
          frames = [frames plant_output_frames.frame];
      else
          frames = [frames plant_output_frames];
      end
      
      checkDependency('lcmgl');
      obj.lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'PFModeHopping');

      obj = setInputFrame(obj,MultiCoordinateFrame.constructFrame(frames));
      obj = setOutputFrame(obj,getStateFrame(obj));
      obj = setSampleTime(obj,getSampleTime(plant));
    end
    
    function x0 = getInitialState(obj)
        x0 = [ones(obj.num_particles, 1)/obj.num_particles; repmat([obj.x0; zeros(obj.nContactForces, 1)], obj.num_particles, 1)];
    end
    
    function inds = getParticleStateInds(obj, i)
        elemsPerState = obj.nx + obj.nContactForces;
        start = obj.num_particles + (1 + elemsPerState*(i-1));
        inds = start:(start + obj.nx - 1);
    end
    function inds = getParticleContactInds(obj, i)
        elemsPerState = obj.nx + obj.nContactForces;
        start = obj.num_particles + (1 + elemsPerState*(i-1)) + obj.nx;
        inds = start:(start + obj.nContactForces - 1);
    end
    
    function xnext = update(obj,t,xin,varargin)
        % if we have inputs parse them now
        % but we don't have inputs yet...
        u = [];
        inputs = obj.getInputFrame.splitCoordinates(varargin{1});
        lidar = reshape(inputs{2}, 3, length(inputs{2})/3);
        kinsol = obj.plant.doKinematics(inputs{1});
        lidar = obj.plant.forwardKin(kinsol, obj.plant.findFrameId('rgbdframe'), lidar);
        
        xnext_pre = xin;
        manip = obj.plant.getManipulator();
        for i=1:obj.num_particles
            inds = obj.getParticleStateInds(i);
            
            %%  PROCESS UPDATE:
            % get basic info about manip in this state
            x = xin(inds);
            q = x(1:obj.nq);
            v = x((obj.nq+1):end);
            c = xin(obj.getParticleContactInds(i));
            [H,C,B] = manipulatorDynamics(manip,q,v);
            
            if (obj.nu > 0)
                tau = B*u-C;
            else
                tau = -C;
            end
            tau = tau + randn(size(tau))*5;
           
            % assemble J
            if obj.nC > 0
                [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q,false);
                
                % Consider updating contact state
                c_proba = exp(-100*(phi+n*v*obj.plant.timestep))*0.9;
                c_norm_active = rand(size(c_proba)) < c_proba;
                
                c_active = zeros(obj.nContactForces, 1);
                c_active(1:(1+obj.nD):end) = c_norm_active;
            
                % construct J and dJ from n,D,dn, and dD so they relate to the
                % lambda vector
                J_norm = zeros(obj.nC, obj.nq);
                J_norm = n;
                
                J_D = zeros(obj.nD*obj.nC, obj.nq);
                for j=1:length(D),
                    J_D(j:obj.nD:end,:) = D{j};
                end
            end
            
            if sum(c_active) > 0
                % solve linear program for contact force, velocity
                % to satisfy (linearized) contact constraints
                % and the dynamics
                
                % decision vars: [qdn; lambda_n]
                
                Aeq = [H/obj.plant.timestep,  -J_norm(c_norm_active~=0, :).';
                     J_norm(c_norm_active~=0, :)*obj.plant.timestep, zeros(sum(c_norm_active~=0), sum(c_norm_active~=0))];
                beq = [tau + H*v/obj.plant.timestep;
                     -phi(c_norm_active~=0)];
                 
                % A*x <= b
                %A = [zeros(sum(c_norm~=0), obj.nq), -eye(sum(c_norm~=0))];
                %b = zeros(sum(c_norm~=0), 1);
                %tmp = quadprog(1E-6*ones(obj.nq+sum(c_norm~=0)),zeros(obj.nq+sum(c_norm~=0), 1), A, b, Aeq, beq);
                
                % for numerical stability:
                Aeq = Aeq + 1E-6*diag(ones(min(size(Aeq)), 1));
                
                tmp = Aeq \ beq;
                vn = tmp(1:obj.plant.getNumPositions);
                qn = q + vn * obj.plant.timestep;
            else
                % no contact forces, so just do standard forward step
                vd = pinv(H)*tau;
                vn = v + vd*obj.plant.timestep;
                qn = q + vn * obj.plant.timestep;
            end
            
            %% MEASUREMENT UPDATE:
            % calculate SDF over the projected location for the lidar
            % points, and sum it up, making weight inverse prop to it.
            kinsol = obj.plant.getManipulator.doKinematics(qn);
            [dists, ~, ~, ~, body_idx] = obj.plant.getManipulator.signedDistances(kinsol, lidar, false);
            sigma = 0.1;
            xnext_pre(i) = sum(exp(-dists(body_idx > 1).^2/sigma^2));
            xnext_pre(inds) = [qn; vn];
            xnext_pre(obj.getParticleContactInds(i)) = c_active;
        end
        
        xnext_pre(1:obj.num_particles)
        
        if (sum(xnext_pre(1:obj.num_particles)) < 1E6)
            xnext_pre(1:obj.num_particles) = xnext_pre(1:obj.num_particles) / sum(xnext_pre(1:obj.num_particles));
        else
            xnext_pre(1:obj.num_particles) = 1 / obj.num_particles;
        end
        
        % Now we need to resample:
        cumu_probs = cumsum(xnext_pre(1:obj.num_particles));
        % trying "stratified resampling" from http://people.isy.liu.se/rt/schon/Publications/HolSG2006.pdf
        cumu_sampling = (rand(obj.num_particles, 1) + (1:obj.num_particles).' - 1)/obj.num_particles;
        xnext = xnext_pre;
        for k=1:obj.num_particles
            ind = find(cumu_probs > cumu_sampling(k), 1);
            xnext(obj.getParticleStateInds(k)) = xnext_pre(obj.getParticleStateInds(ind));
            xnext(obj.getParticleContactInds(k)) = xnext_pre(obj.getParticleContactInds(ind));
        end
        
        if ~isa(obj.v, 'PlanarRigidBodyVisualizer')
            kinsol = obj.plant.doKinematics( inputs{1} );
            for i = 2:(length(obj.plant.getManipulator.body))
                bdy = obj.plant.getManipulator.body(i);
                obj.lcmgl.glColor4f(0, 1, 1, 0.1);
                posrpy = obj.plant.forwardKin(kinsol, i, bdy.visual_geometry{1}.T(1:3, 4), 1);
                axis_angle = rpy2axis(posrpy(4:6));
                obj.lcmgl.glPushMatrix();
                obj.lcmgl.glTranslated(posrpy(1), posrpy(2), posrpy(3));
                obj.lcmgl.glRotated(axis_angle(4)*180/pi, axis_angle(1), axis_angle(2), axis_angle(3));
                obj.lcmgl.box(zeros(3, 1), bdy.visual_geometry{1}.size);
                obj.lcmgl.glPopMatrix();
            end
        end

    end
    
    function y = output(obj,t,x,varargin)
      if (t > 0.01)
          inputs = obj.getInputFrame.splitCoordinates(varargin{1});
          lidar = reshape(inputs{2}, 3, length(inputs{2})/3);
          kinsol = obj.plant.doKinematics(inputs{1});
          lidar = obj.plant.forwardKin(kinsol, obj.plant.findFrameId('rgbdframe'), lidar);
          
          if isa(obj.v, 'PlanarRigidBodyVisualizer')
            clf; hold on;
          end
          % why oh why is u (varargin) all zeros.
          % whyyyy
          %scatter(lidar(1, :), lidar(3, :));
          
          for k = 1:obj.num_particles
            inds = obj.getParticleStateInds(k);
            if isa(obj.v, 'PlanarRigidBodyVisualizer')
                obj.v.drawBody(t, x(inds(1:obj.nq)), 1, min(x(k)*5, 1.0), [1, 0, 1]);
            else
                kinsol = obj.plant.doKinematics( x(obj.getParticleStateInds(k)) );
                for i = 2:(length(obj.plant.getManipulator.body))
                    bdy = obj.plant.getManipulator.body(i);
                    obj.lcmgl.glColor4f(1, 0, 1, min(x(k)*1, 1.0));
                    posrpy = obj.plant.forwardKin(kinsol, i, bdy.visual_geometry{1}.T(1:3, 4), 1);
                    axis_angle = rpy2axis(posrpy(4:6));
                    obj.lcmgl.glPushMatrix();
                    obj.lcmgl.glTranslated(posrpy(1), posrpy(2), posrpy(3));
                    obj.lcmgl.glRotated(axis_angle(4)*180/pi, axis_angle(1), axis_angle(2), axis_angle(3));
                    obj.lcmgl.box(zeros(3, 1), bdy.visual_geometry{1}.size);
                    obj.lcmgl.glPopMatrix();
                end
            end
          end

          obj.lcmgl.switchBuffers();
          
      end
      
      %y = zeros(obj.nx, 1);
      %for i=1:obj.num_particles
      %   y = y + x(i) * x(obj.getParticleStateInds(i)); 
      %end
      %x(1:obj.num_particles);
      y = x;
      
    end
  end
  
end