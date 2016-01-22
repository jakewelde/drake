classdef PFModeHopping < DrakeSystem
    
  properties
      x0; % initial state of the system we're estimating
      nx; % state size of system we're estimating
      nu; % inputs of system we're estimating
      nz; % outputs of system we're estimating
      plant; % system we're estimating
      v; % viewer of it
      
      num_particles;
      
      % STATE VECTOR:
      % [ particle_weights; % num_particles x 1, should be normalized
      %   particle_1_state; % nX x 1
      %   ...
      %   particle_<num_particles>_state; % nX x 1
      % ]
  end
  
  methods
    function obj = PFModeHopping(plant, num_particles, x0, state_frame)
      typecheck(plant,'TimeSteppingRigidBodyManipulator');
      
      nx = getNumStates(plant);
      nu = getNumInputs(plant);
      nz = getNumOutputs(plant);
      num_states = num_particles + num_particles*nx;
      obj = obj@DrakeSystem(0,num_states,nu+nz,nx,0,1);
      obj.num_particles = num_particles;
      
      obj.x0 = x0;
      obj.plant = plant;
      obj.v = obj.plant.constructVisualizer();
      obj.nx = nx;
      obj.nu = nu;
      obj.nz = nz;
      
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
        x0 = [ones(obj.num_particles, 1)/obj.num_particles; repmat(obj.x0, obj.num_particles, 1)];
    end
    
    function inds = getParticleStateInds(obj, i)
        inds = obj.num_particles+1 + obj.nx*(i-1): obj.num_particles + obj.nx*i;
    end
    
    function xnext = update(obj,t,x,varargin)
        % if we have inputs parse them now
        % but we don't have inputs yet...
        xnext = x;
        for i=1:obj.num_particles
           inds = obj.getParticleStateInds(i);
           xnext(inds) = obj.plant.update(t, x(inds), []) + randn(obj.nx, 1)*0.01;
        end
    end
    
    function y = output(obj,t,x,varargin)
      if (t > 0.01)
          clf; hold on;
          for k = 1:obj.num_particles
            inds = obj.getParticleStateInds(k);
            obj.v.drawBody(t, x(inds(1:obj.plant.getNumPositions)), 1);
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