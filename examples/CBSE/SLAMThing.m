classdef SLAMThing < MIMODrakeSystem
  
  properties
    t_history;
    state_history;
    imu_history;
    num_potential_contacts;
    r;
    dt;
    imu;
    last_sol;
    last_R;
    last_d;
    
    % for ze drawing
    lcmgl;
    
    % contact info
    nC; % # contact pairs
    nD; % 
    active_collision_options;
    nContactForces; 
    
    eval_time;
    
  end
  
  methods
    function obj = SLAMThing(r,options)
      if nargin<2
        options = struct();
      end
      
      input_frame = getOutputFrame(r);
      output_frame = getOutputFrame(r);
      
      output_frame = {r.getStateFrame, r.getOutputFrame.getFrameByName('brickIMU'), r.getStateFrame};
      output_frame = MultiCoordinateFrame(output_frame);
      
      obj = obj@MIMODrakeSystem(0,0,input_frame,output_frame,true,true);
      obj = setInputFrame(obj,input_frame);

      %obj = obj.setNumOutputs(length(output_frame.getCoordinateNames));
      obj = setOutputFrame(obj,output_frame);
      
      % should take from options
      dt = 0.01;
      obj.dt = dt;
      obj = setSampleTime(obj,[0.01;0]); % sets controller update rate
      
      % Setup memory of state and imu
      obj.t_history = SharedDataHandle([]);
      obj.state_history = SharedDataHandle([]);
      obj.imu_history = SharedDataHandle([]);
      obj.last_sol = SharedDataHandle([]);
      obj.last_R = SharedDataHandle([]);
      obj.last_d = SharedDataHandle([]);
      
      obj.eval_time = SharedDataHandle([]);
      
      % Store robot
      obj.r = r;
      % Currently ONLY thinking about terrain... need to figure out the
      % way contacts are dealt with in general.
      %obj.num_potential_contacts = size(r.getTerrainContactPoints, 2);
      
      % Generate an imu for us to use as measurement func
      body = 2;
      xyz = [0; 0; 0];
      rpy = zeros(3, 1); %nonzero rpy not supported in imu
      obj.imu = RigidBodyInertialMeasurementUnit(r, body, xyz, rpy);
      
      % be able to draw stuff!
      checkDependency('lcmgl');
      obj.lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'SLAM!');
      
      % set up some utility stuff to help keep track of contacts
      obj.nC = obj.r.getManipulator.getNumContactPairs;
      [~,normal,d] = obj.r.getManipulator.contactConstraints(zeros(obj.r.getManipulator.getNumPositions,1));
      obj.nD = 2*length(d);
      assert(size(normal,2) == obj.nC); % just a double check    
      obj.active_collision_options.terrain_only = true;
      obj.nContactForces = obj.nC*(2 + obj.nD); 
      
    end
    
    function [utility, j] = objective_fun(obj, x)
      % Given a state, predicts what the imu output will be.
      imu_hist = obj.imu_history.getData;
      imu_pred = imu_hist * 0;
      j = sparse(x.'*0);
      for i=1:length(x)/12
        % Transform to local frame to get im
        
        % Predicted vals for omega
        offset = 12*(i-1);
        i_range = 1+offset:12+offset;
        x_i = x(i_range);
        
        %imu_pred(:, i) = obj.imu.output(obj.r.getManipulator,0,double(x_i),[]);
        [out, acc, gyro] = obj.imu.output_and_jacobian(obj.r.getManipulator,0,double(x_i), []);
        % Numerical:
        %         out_norm = obj.imu.output(obj.r.getManipulator,0,double(x_i), []);
        %         out_num = zeros(6, length(x_i));
        %         eps = 1E-3;
        %         for j=1:length(x_i)
        %           add = x_i*0;
        %           add(j) = eps;
        %           out_this = (obj.imu.output(obj.r.getManipulator,0,double(x_i+add), []) - out_norm)/eps;
        %           acc = out_this(8:10);
        %           gyro = out_this(5:7);
        %           out_num(:, j) = [acc; gyro];
        %         end
        %         acc = out_norm(8:10);
        %         gyro = out_norm(5:7);
        %         out = out_num;
        
        
        imu_pred(5:7, i) = gyro;
        imu_pred(8:10, i) = acc;
        error = imu_hist(5:10, i) - imu_pred(5:10, i);
        j(i_range) = out.'*error / obj.dt;
      end
      % Total error, then, is...
      err = imu_pred(5:10, :) - imu_hist(5:10, :);
      j = [];
      utility = sum(diag(err.' * err));
    end
    
    function [integration_constraint, integration_equality, manip_constraint, manip_equality, imu_constraint, imu_equality] = ISAM_generate_constraints(obj, n, nX, m)
      state_hist = obj.state_history.getData;
      imu_hist = obj.imu_history.getData;
      % (6 dynamic eqs, 6 manip eqs, 6 measurement)
      % (when n=1, only measurements present)
      
      nX_per = 12 + m;
 
      if (n~=1)
        % Generate dynamic constraints
        % velocity integrates into position for each time point
        integration_constraint = zeros(6, nX);
        integration_equality = zeros(6, 1);
        offset = nX_per*(n-2);
        integration_constraint(1:6, 1+offset:6+offset) = eye(6); %q_n
        integration_constraint(1:6, nX_per+1+offset:nX_per+6+offset) = -eye(6); %-q_n+1
        integration_constraint(1:6, nX_per+7+offset:nX_per+12+offset) = obj.dt*eye(6);%+h*qd_n+1 (backwards integration here)

        % And the manipulator equation
        manip_constraint = zeros(6, nX);
        manip_equality = zeros(6, 1);
        [Hnext,Cnext,Bnext] = obj.r.manipulatorDynamics(state_hist(1:6, n), state_hist(7:end, n));
        
        [phi,normal,~,~,~,~,~,~,nc,D,dnc,dD] = obj.r.contactConstraints(state_hist(1:6, end),false, obj.active_collision_options);
        % construct J and dJ from nc,D,dnc, and dD so they relate to the
        % lambda vector
%         J = zeros(m,6);
%         J(1:2+obj.nD:end,:) = n;
%         dJ = zeros(m*6,6);
%         dJ(1:2+obj.nD:end,:) = dnc;
%         for j=1:length(D),
%           J(1+j:2+obj.nD:end,:) = D{j};
%           dJ(1+j:2+obj.nD:end,:) = dD{j};
%         end
          
        %[~,Gnext] = obj.r.manipulatorDynamics(state_hist(1:6, n),zeros(6, 1));
        offset = nX_per*(n-2);
        % Hnext * (qd_n+1 - qd) +
        % dt*(C_next+Gnext-B_next*u_next-J_next.'*contact_forces)
        % So dropping some terms for now (no u, no contact forces)
        manip_constraint(:, nX_per+7+offset:nX_per+12+offset) = Hnext;
        manip_constraint(:, 7+offset:12+offset) = -Hnext;
        % -J.' * contact forces
        %manip_constraint(:, nX_per+13+offset:nX_per+12+m+offset) = -J.';
        
        manip_equality = manip_equality - obj.dt * Cnext;%; - obj.dt*Gnext;
        
      else
        integration_constraint = [];
        integration_equality = [];
        manip_constraint = [];
        manip_equality = [];
      end
      
      % And measurement constraint on IMU
      % -> given J_imu at x_n, J_imu(x_n) * (x_n+1 - x_n) + J0 = z_n+1
      imu_constraint = zeros(6, nX);
      imu_equality = zeros(6, 1);
      x = obj.last_sol.getData;
      if (isempty(x))
        x = zeros(12+m, 1);
      else
        x = [x; x(end-m-11:end-m)];
      end
      offset = nX_per*(n-1);
      i_range = 1+offset:12+offset;
      x_i = x(i_range);
      % Analytic (probably not right D:):
      %out = obj.imu.output_and_jacobian(obj.r.getManipulator,0,double(x_i), []);
      % Numerical:
      out_norm = obj.imu.output(obj.r.getManipulator,0,double(x_i), []);
      out_num = zeros(6, length(x_i));
      eps = 1E-6;
      add = x_i*0;
      for j=1:length(x_i)
        add = add*0;
        add(j) = eps;
        out_this = (obj.imu.output(obj.r.getManipulator,0,double(x_i+add), []) - out_norm)/eps;
        acc = out_this(8:10);
        gyro = out_this(5:7);
        out_num(:, j) = [acc; gyro];
      end
      acc = out_norm(8:10);
      gyro = out_norm(5:7);
      out = out_num;

      % 1st-order taylor expansion of imu output at x_n.
      % J * (xn+1 - x_n) + zn = zn+1
      % i.e. j*xn+1 = zn+1 - J0 + J*x_n
      % Jac * x_n+1 - Jac * x_n
      
      if (n~=1)
        imu_constraint(:, offset-m-11:offset-m) = -out;
        imu_constraint(:, offset+1:offset+12) = out;
      end
      %imu_constraint(offset_inner+1:offset_inner+6, offset+1:offset+12) = -out;
      % ... J0 = z_n+1
      %imu_constraint(offset_inner+1:offset_inner+6, offset+13:offset+18) = eye(6);
      imu_equality(1:3) = imu_hist(8:10, n) - acc;
      imu_equality(4:6) = imu_hist(5:7, n) - gyro;
      %imu_equality = imu_equality; % + out*double(x_i);
    end
    
    function doISAM(obj)
      state_hist = obj.state_history.getData;
      % Decision variables: given N time steps so far,
      % N*12 (x and x* states) %% (later) + m*N (m contact forces)
      N = size(state_hist, 2);
      m = 0; %obj.nContactForces;
      nX = N*12 + m*N;
      % Constraints:
      % (N-1)*12 + N*6 (6 dynamic eqs, 6 manip eqs, 6 measurement)
      nCons = (N-1)*12 + N*6;
      A = zeros(nCons, nX);
      b = zeros(nCons, 1);
      
      [integration_constraint, integration_equality, ...
        manip_constraint, manip_equality, ...
        imu_constraint, imu_equality] = ISAM_generate_constraints(obj, N, nX, m);
      if (N~=1)
        A_next = [integration_constraint; manip_constraint; imu_constraint];
        b_next = [integration_equality; manip_equality; imu_equality];
      else
        A_next = imu_constraint;
        b_next = imu_equality;
      end
      A(end-size(A_next, 1)+1:end, :) = A_next;
      b(end-size(b_next, 1)+1:end) = b_next;
      
      % Solve
      if size(A, 1) > 0
        
        % If we don't currently have a R and d established,
        % or we want to regenerate (and re-linearize),
        % generate them the long way.
        R = sparse(obj.last_R.getData());
        d = sparse(obj.last_d.getData());
        if (isempty(d) || isempty(R) || mod(N, 25) == 0)
          % We need to fill in the rest of A, then...
          for n=1:N-1
            [integration_constraint, integration_equality, ...
              manip_constraint, manip_equality, ...
              imu_constraint, imu_equality] = ISAM_generate_constraints(obj, n, nX, m);
            if (n~=1)
              A_next = [integration_constraint; manip_constraint; imu_constraint];
              b_next = [integration_equality; manip_equality; imu_equality];
              A(1 + (n-2)*12+(n-1)*6 : (n-1)*12+(n)*6, :) = A_next;
              b(1 + (n-2)*12+(n-1)*6 : (n-1)*12+(n)*6) = b_next;
            else
              A_next = imu_constraint;
              b_next = imu_equality;
              A(1:6, :) = A_next;
              b(1:6) = b_next;
            end
          end
          % For variables that are fully undetermined (e.g. yaw),
          % add some regularization (force to zero).
          %Areg = [A; 1E-6 * eye(size(A, 2))];
          %breg = [b; zeros(size(A, 2), 1)];
          
          % Bind the first position to the origin
          Areg = [zeros(12+m, size(A, 2)); A];
          Areg(1:12+m, 1:12+m) = eye(12+m);
          binit = [ 0
                    0
                    10
                    0.5
                    0
                    0
                    0
                    0
                    0
                    0
                    0
                    0
                    zeros(m, 1)];
          breg = [binit; b];
          
          [Q,R] = qr(Areg);
          R = R(1:size(R, 2), :);
          de = Q.'*breg;
          d = de(1:size(R, 1));
        else
          %Otherwise qr-update!
          new_A = [manip_constraint(end-5:end, :);
            integration_constraint(end-5:end, :);
            imu_constraint(end-5:end, :)];
          new_b = [manip_equality(end-5:end);
            integration_equality(end-5:end);
            imu_equality(end-5:end)];
          
          % Don't forget
          newR = sparse(zeros(size(R, 1)+size(new_A, 1), size(new_A, 2)));
          newR(1:size(R, 1), 1:size(R,2)) = R;
          R = newR;
          R(end-size(new_A, 1)+1:end, 1:size(new_A, 2)) = new_A;
          d = [d; new_b];
          
          % Proceeding columnwise through all zeros...
          [i, k] = find(R);
          lower_triangle = k < i;
          i = i(lower_triangle);
          k = k(lower_triangle);
          upper = length(i);
          for j=1:upper
            % only those below diagonal
            [c s] = givens_isam(R(k(j), k(j)), R(i(j), k(j)));
            %mod_mat = sparse(eye(size(R, 1)));
            % Don't know why, but this doesn't work if this isn't
            % negative of the formula described in Kaess '08.
            %               mod_mat(i(j), i(j)) = -c;
            %               mod_mat(k(j), k(j)) = -c;
            %               mod_mat(i(j), k(j)) = -s;
            %               mod_mat(k(j), i(j)) = s;
            
            new_row_k = -c*R(k(j), :) + -s*R(i(j), :);
            new_row_i = s*R(k(j), :) + -c*R(i(j), :);
            new_row_k(abs(new_row_k)<1E-5) = 0;
            new_row_i(abs(new_row_i)<1E-5) = 0;
            R(k(j), :) = new_row_k;
            R(i(j), :) = new_row_i;
            
            new_row_k = -c*d(k(j)) + -s*d(i(j));
            new_row_i = s*d(k(j)) + -c*d(i(j));
            d(k(j)) = new_row_k;
            d(i(j)) = new_row_i;
            
            %R = mod_mat*R;
            %d = mod_mat*d;
          end
          %R(abs(R)<1E-5) = 0;
        end
        
        theta = R\d;
        %           elems = size(R, 2);
        %           theta = zeros(elems, 1);
        %           for ii = elems : -1 : 1
        %               theta(ii,:) = (d(ii,:)-R(ii,:)*theta)/R(ii,ii);
        %           end
        
        obj.last_sol.setData(theta);
        obj.last_R.setData(R);
        obj.last_d.setData(d);
                    q_guess = [theta(end-m-11:end-m-6)]
                    q_d_guess = theta(end-m-5:end-m)
        error = theta(end-m-11:end-m) - state_hist(:, end)
        N
        if ~isempty(obj.lcmgl)
          obj.lcmgl.glColor3f(1, 0, 0);

          obj.lcmgl.glTranslated(theta(end-11), theta(end-10), theta(end-9)); 
          
          obj.lcmgl.glRotated(theta(end-8)*180/pi, 1, 0, 0);
          obj.lcmgl.glRotated(theta(end-7)*180/pi, 0, 1, 0);
          obj.lcmgl.glRotated(theta(end-6)*180/pi, 0, 0, 1);
          obj.lcmgl.box([0,0,0], [2, 0.5, 1]);
          
          
  %          obj.lcmgl.line3(origin(1), origin(2), origin(3), point(1,i), point(2,i), point(3,i));
          obj.lcmgl.switchBuffers;
        end
      
      else
        obj.last_sol.setData(zeros(nX, 1));
      end
    end
    
    function varargout=mimoOutput(obj,t,~,varargin)
      tic
      t_hist = obj.t_history.getData;
      % some logic to skip the first time we're called when t=0,
      % but remember the second time (when the state is properly set up)
      repeat = false;
      if (t == 0)
        if (t == obj.t_history.getData())
          repeat = true;
        end
        obj.t_history.setData(t);
      end
      % only update based on our own dt, which might subsample
      % true update rate... hack, should do cleaner, I'm sure drake
      % can do this internally
      if ((t == 0 && repeat) || (t > 0 && (isempty(t_hist) || (t - t_hist(end) >= obj.dt-1E-5))))
        % Store this state persistently
        obj.t_history.setData([obj.t_history.getData t]);
        obj.state_history.setData([obj.state_history.getData varargin{1}]);
        obj.imu_history.setData([obj.imu_history.getData varargin{2}]);
        %obj.doNonlinOptim();
        obj.doISAM();
        
      end
      last_sol = obj.last_sol.getData();
      if (isempty(last_sol))
        last_sol = zeros(length(obj.getOutputFrame.getFrameByNum(3).getCoordinateNames), 1);
      else
        last_sol = full(last_sol(end-11:end));
      end
      varargout = [varargin last_sol];
      tim = toc
      
      obj.eval_time.setData([obj.eval_time.getData; tim]);
      obj.eval_time.getData
    end
  end
  
end
