classdef QPLCContactImplicitEstimator < DrakeSystem
    
  properties
      x0;
      nu;
      nx;
      nz;
      plant;
      v;
      sensor_data;
      
      lcmgl;
      
      
      use_nlp;
  end
  
  methods
    function obj = QPLCContactImplicitEstimator(plant, x0, sensor_data)
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
      u = w(1:obj.nu);
      x_gt = w(obj.nu+1:obj.nu+obj.nx);
      z = reshape(w(obj.nu+obj.nx+1:end), 3, (length(w)-obj.nx-obj.nu)/3);
      dt = obj.getSampleTime();
      dt = dt(1);
      
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
      dt = dt / iters;
      xlast = x;
      for iter=1:iters
          qlast = xlast(1:nq);
          vlast = xlast(nq+1:end);
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
        
      
          % Solve for positions that explain the sensor data as a baseline
          if (obj.use_nlp)
              nVars = nq;
              solver = NonlinearProgram(nVars);
              sdfCost = FunctionHandleConstraint(0,0,nq, @(qhat)obj.sdfObjective(qhat,z));
              solver = solver.addCost(sdfCost, 1:nVars);
              solver = solver.setSolverOptions('snopt', 'MajorOptimalityTolerance', 1E-3);
              solver = solver.setSolver('fmincon');
              solver = solver.addDisplayFunction(@(qhat)obj.v.drawWrapper(0, qhat));
              [qlast, F, info, infeasible_constraint_name] = solver.solve(qlast);
              xlast = [qlast; zeros(nq, 1)]
          else
              % as a local approximation, project all cloud points onto the
              % surface of the object positioned via the LAST state estimate
             
              % can we vectorize better?    
              f = zeros(nq, 1);
              Q = zeros(nq, nq);
              K = 0;
              unique_bodies = unique(body_idx);
              pointcloud_weight = 1.0 / numel(body_idx);
              for k=1:numel(unique_bodies)
                  this_body = body_idx==unique_bodies(k);
                  on_this_body = sum(this_body);
                  [~, J_zi] = obj.plant.forwardKin(kinsol, unique_bodies(k), body_z_prime(:, this_body));
                  Ks = z(:, this_body) - z_prime(:, this_body) + reshape(J_zi*qlast, 3, on_this_body);
                  f = f - pointcloud_weight*(2 * reshape(Ks, on_this_body*3, 1).' * J_zi).';
                  Q = Q + pointcloud_weight*2*(J_zi.' * J_zi).';
                  K = K + pointcloud_weight*sum(sum(Ks.*Ks));
              end
              
              euler_weight = 0.1;
              target_q = qlast + dt * vlast;
              f = f - euler_weight * (2 * target_q);
              Q = Q + euler_weight * (2 * eye(nq));
              K = K + euler_weight * (target_q.' * target_q);
              
              % If Q isn't full rank we'll have to hack it to be...
              while min(eig(Q)) < 0
                  'insufficient data'
                  [V,D] = eig(Q);
                  Q = Q + V(:, 1)*V(:, 1).'*(eps(D(1,1))-D(1,1));
              end

    %           % or, not vectorized:          
    %           f = zeros(6, 1);
    %           Q = zeros(6, 6);
    %           for i=1:size(z, 2)
    %               [~, J_zi] = obj.plant.forwardKin(kinsol, body_idx(i), body_z_prime(:, i));
    %               K = z(:, i) - z_prime(:, i) + J_zi * qlast;
    %               f = f - (2 * K.'*J_zi).';
    %               Q = Q + 2*(J_zi.' * J_zi).';
    %           end

              solver = QuadraticProgram(Q, f);
              [qlast, F, info, active] = solver.solve(qlast);
              F = F + K
              xlast = [qlast; x_gt(nq+1:end)]
              
              obj.v.draw(0, xlast(1:nq));
          end
          xnext = xlast;
      end
      
    end
    
    function [c, dc] = sdfObjective(obj, q, z)
        kinsol = obj.plant.doKinematics(q);
        [phi,n,x,x_body,body_idx] = obj.plant.getManipulator.signedDistances(kinsol,z,false);
        K = 1000 / size(z, 2);
        c = K*phi.'*phi;
        
        % derivatives
        bodies = unique(body_idx);
        dc = zeros(1,obj.plant.getNumPositions);
        for body=bodies.'
            [~, J] = obj.plant.forwardKin(kinsol, body, x_body);
            for i = find(body_idx == body).'
                dc = dc - K * 2 * phi(i) * n(:,i).'*J((i-1)*3+1:i*3, :);
            end
        end
    end
    
    function y = output(obj,t,x,varargin)
      y = x;
    end
    
  end
  
end