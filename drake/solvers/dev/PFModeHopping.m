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
      obj = obj@DrakeSystem(0,num_states,nu+nz,nx,0,1);
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

      obj = setInputFrame(obj,MultiCoordinateFrame.constructFrame(frames));
      obj = setOutputFrame(obj,getStateFrame(plant));
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
        xnext = xin;
        manip = obj.plant.getManipulator();
        for i=1:obj.num_particles
            inds = obj.getParticleStateInds(i);
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
            
            % Consider updating contact state
            c = rand(size(c)) > 0.5;
           
            % assemble J
            if obj.nC > 0
                [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q,false);
                % construct J and dJ from n,D,dn, and dD so they relate to the
                % lambda vector
                J = zeros(obj.nContactForces,obj.nq);
                J(1:1+obj.nD:end,:) = n;
                %dJ = zeros(obj.nContactForces*nq,nq);
                %dJ(1:1+obj.nD:end,:) = dn;
                
                for j=1:length(D),
                    J(1+j:1+obj.nD:end,:) = D{j};
                    %dJ(1+j:1+obj.nD:end,:) = dD{j};
                end
                
                %fv = fv + J'*lambda;
                %dfv(:,1:nq) = dfv(:,1:nq) + matGradMult(dJ,lambda,true);
                %dfv(:,1+nq+nu:nq+nu+nl) = J'*obj.options.lambda_mult;
            end
            
            activeJ = J(c~=0, :);
            
            if min(size(activeJ)) > 0
                % do an optimization to produce contact forces to avoid
                % penetration
                contactForces = rand(sum(c), 1);
                tau = tau + activeJ.' * contactForces;
            end
      
           vd = pinv(H) * tau + pinv(H) * tau .* randn(size(v))*0.05;
           q = q + v * obj.plant.timestep;
           v = v + vd * obj.plant.timestep;
           
           xnext(inds) = [q; v];
        end
    end
    
    function y = output(obj,t,x,varargin)
      if (t > 0.01)
          clf; hold on;
          for k = 1:obj.num_particles
            inds = obj.getParticleStateInds(k);
            obj.v.drawBody(t, x(inds(1:obj.nq)), 1);
          end
      end
      
      y = zeros(obj.nx, 1);
      for i=1:obj.num_particles
         y = y + x(i) * x(obj.getParticleStateInds(i)); 
      end
      y
      
    end
  end
  
end