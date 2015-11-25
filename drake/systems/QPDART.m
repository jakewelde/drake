classdef QPDART < DrakeSystem
    
  properties
      x0;
      nu;
      nx;
      nz;
      nC;
      nD;
      plant;
      v;
      sensor_data;
      
      lcmgl;
      
      
      use_nlp;
      
      point_to_plane = 0;
  end
  
  methods
    function obj = QPDART(plant, x0, sensor_data)
      typecheck(plant,'TimeSteppingRigidBodyManipulator');
      
      obj = obj@DrakeSystem(0,plant.getNumStates(),plant.getNumInputs()+plant.getNumOutputs(),plant.getNumStates(),0,1);

      obj.x0 = x0;
      obj.sensor_data = sensor_data;
      
      obj.use_nlp = 0;
      
      % setup input/output frames
      plant_input_frames = getInputFrame(plant);
      obj.nu = plant.getNumInputs;
      plant_output_frames = getOutputFrame(plant);
      obj.nz = plant.getNumOutputs;
      obj.nx = plant.getNumStates;
      obj.nC = plant.getNumContactPairs;
      [~,~,d] = plant.contactConstraints(zeros(plant.getNumPositions,1));
      obj.nD = 2*length(d);
      obj.plant = plant;
      
      obj.v = obj.plant.constructVisualizer();
      
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
      
      checkDependency('lcmgl');
      obj.lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'case');
      
    end
    
    function x0 = getInitialState(obj)
        x0 = obj.x0;
    end
    
    function xnext = update(obj,t,x,varargin)
      w = varargin{1};
      x_gt = w(obj.nu+1:obj.nu+obj.nx);
      z = reshape(w(obj.nu+obj.nx+1:end), 3, (length(w)-obj.nx-obj.nu)/3);
      
      nq = obj.plant.getNumPositions;
      
      
      % z is now raycast points in world frame
      ds = obj.plant.getManipulator.sensor{2};
      kinsol = obj.plant.doKinematics(x(1:nq));
      z = obj.plant.forwardKin(kinsol, ds.frame_id, z);
      % throw out ones that are too far away from ground truth
      kinsol_gt = obj.plant.doKinematics(x_gt(1:nq));
      [phi,~,~,~,body_idx] = obj.plant.getManipulator.signedDistances(kinsol_gt,z,false);
      
      keep = phi < 0.1 & body_idx ~= 1;
      z = z(:, keep);
      
      iters = 1;
      xlast = x;
      zlast = [xlast(1:nq)];
      for iter=1:iters
          qlast = xlast(1:nq);
          kinsol = obj.plant.doKinematics(qlast);
          [phi,n,z_prime,body_z_prime,body_idx] = obj.plant.getManipulator.signedDistances(kinsol,z,false);
          
          
          obj.lcmgl.glColor3f(1, 0, 0);
          obj.lcmgl.points(z(1,:),z(2,:),z(3,:));
          obj.lcmgl.glColor3f(0, 1, 0);
          obj.lcmgl.points(z_prime(1,:),z_prime(2,:),z_prime(3,:));
          for i=1:size(z, 2)
             obj.lcmgl.line3(z(1,i),z(2,i), z(3,i), z_prime(1,i),z_prime(2,i),z_prime(3,i)); 
          end
          obj.lcmgl.switchBuffers;
        
      
          % as a local approximation, project all cloud points onto the
          % surface of the object positioned via the LAST state estimate

          % can we vectorize better?    
          nvars = nq; % + obj.nC*(obj.nD+1);
          qinds = 1:nq;

          f = zeros(nvars, 1);
          Q = zeros(nvars, nvars);
          K = 0;

          unique_bodies = unique(body_idx);
          pointcloud_weight = 1.0 / numel(body_idx);
          for k=1:numel(unique_bodies)
              this_body = body_idx==unique_bodies(k);
              this_body_i = find(this_body);
              on_this_body = sum(this_body);
              [~, J_zi] = obj.plant.forwardKin(kinsol, unique_bodies(k), body_z_prime(:, this_body));
              Ks = z(:, this_body) - z_prime(:, this_body) + reshape(J_zi*qlast, 3, on_this_body);
              if (obj.point_to_plane)
                  for l=1:on_this_body
                      ns = n(:, this_body_i(l));
                      J_zil = J_zi((3*(l-1)+1):3*l, :);
                      f(qinds) = f(qinds) - pointcloud_weight*(2 * Ks(:, l).' * (ns * ns.') * J_zil).';
                      Q(qinds, qinds) = Q(qinds, qinds) + pointcloud_weight*2*(J_zil.' * (ns * ns.') * J_zil).';
                      K = K + pointcloud_weight*Ks(:, l).'*Ks(:, l);
                  end
              else
                  f(qinds) = f(qinds) - pointcloud_weight*(2 * reshape(Ks, on_this_body*3, 1).' * J_zi).';
                  Q(qinds, qinds) = Q(qinds, qinds) + pointcloud_weight*2*(J_zi.' * J_zi).';
                  K = K + pointcloud_weight*sum(sum(Ks.*Ks));
              end
          end
          
          % If Q isn't full rank we'll have to hack it to be...
          % fill in unspecified cols of Q
          Q = Q + diag(sum(Q)==0);
          i = 1000;
          while min(eig(Q)) < 0 && i > 0
              if i == 1000
                'insufficient data'
              end
              [V,D] = eig(Q);
              Q = Q + V(:, 1)*V(:, 1).'*(eps(D(1,1))-D(1,1));
              i = i - 1;
          end
          
          solver = QuadraticProgram(Q, f);
          [zlast, F, info, active] = solver.solve(zlast);
          qlast = zlast(qinds);
          xlast = [qlast; 0*qlast]
          F = F + K

          obj.v.draw(0, qlast);
          
          xnext = xlast;
      end
      
    end
    
    function y = output(obj,t,x,varargin)
      y = x;
    end
    
  end
  
end