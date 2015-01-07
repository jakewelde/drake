classdef StateEstQPBlock < MIMODrakeSystem
  
  properties
    t_history;
    state_history;
    imu_history;
    r;
    dt;
    imu;
    last_sol;
    full_state_est_tt;
    
    % subframe of input state frame representing our robot itself
    state_frame_subind;
    
    % for ze drawing
    lcmgl;
    
    eval_time;
    init_cond;
    
    % Cache various constraint matrices
    F_manip_cache;
    g_manip_cache;
    F_meas_cache;
    g_meas_cache;
    A_int_cache;
    b_int_cache;
    sense_int_cache;
    A_nonpenetration_cache;
    b_nonpenetration_cache;
    sense_nonpenetration_cache;
      
    A_init; 
    b_init;
    sense_init;
    
  end
  
  methods
    function obj = StateEstQPBlock(r,init_cond, options, ...
      full_state_est_tt_handle, full_state_est_handle)
      if nargin<2
        init_cond = zeros(12, 1);
      end
      if (nargin<3)
        options = struct();
      end
      if (nargin<4)
        full_state_est_tt_handle = SharedDataHandle([]);
      end
     if (nargin<5)
        full_state_est_handle = SharedDataHandle([]);
      end
      if ~isfield(options,'dt_state_est')
        options.dt_state_est = 0.01;
      end
      
      if (~isfield(options,'state_frame_subind'))
        options.state_frame_subind = 0;
        output_frame = {r.getStateFrame, r.getOutputFrame.getFrameByName('brickIMU'), r.getStateFrame};
      else
        robot_state = r.getStateFrame.getFrameByNum(options.state_frame_subind);
        output_frame = {robot_state, r.getOutputFrame.getFrameByName('brickIMU'), robot_state};
      end
      
      input_frame = getOutputFrame(r);
      output_frame = MultiCoordinateFrame(output_frame);
      
      obj = obj@MIMODrakeSystem(0,0,input_frame,output_frame,true,true);
      
      % take things out of options
      obj.dt = options.dt_state_est;
      obj = setSampleTime(obj,[obj.dt;0]); % sets controller update rate
      obj.state_frame_subind = options.state_frame_subind;
      obj.init_cond = init_cond;
      
      % Can set up initial condition matrices easily
      % Only supply init conditions for pos, linear velocity, and yaw
      obj.A_init = sparse(12, 12); 
      obj.b_init = zeros(12, 1);
      obj.sense_init = zeros(12, 1);
%       obj.A_init(1:3, 1:3) = eye(3);
%       obj.b_init(1:3) = obj.init_cond(1:3);
%       obj.A_init(6, 6) = 1;
%       obj.b_init(6) = obj.init_cond(6);
%       obj.A_init(7:9, 7:9) = eye(3);
%       obj.b_init(7:9) = obj.init_cond(7:9);
      obj.A_init = eye(12);
      obj.b_init = obj.init_cond;
      obj.sense_init(:) = '=';
      
      % Setup memory of state and imu
      obj.t_history = SharedDataHandle([]);
      obj.state_history = SharedDataHandle([]);
      obj.imu_history = SharedDataHandle([]);
      obj.last_sol = full_state_est_handle;
      obj.full_state_est_tt = full_state_est_tt_handle;
      obj.eval_time = SharedDataHandle([]);
      
      % And caches for various matrices
      obj.F_manip_cache = SharedDataHandle([]);
      obj.g_manip_cache = SharedDataHandle([]);
      obj.F_meas_cache = SharedDataHandle([]);
      obj.g_meas_cache = SharedDataHandle([]);
      obj.A_int_cache = SharedDataHandle([]);
      obj.b_int_cache = SharedDataHandle([]);
      obj.sense_int_cache = SharedDataHandle([]);
      obj.A_nonpenetration_cache = SharedDataHandle([]);
      obj.b_nonpenetration_cache = SharedDataHandle([]);
      obj.sense_nonpenetration_cache = SharedDataHandle([]);

      % Store robot
      obj.r = r;
      
      % Generate an imu for us to use as measurement func
      body = 2;
      xyz = [0; 0; 0];
      rpy = zeros(3, 1); %nonzero rpy not supported in imu
      obj.imu = TimeSteppingInertialMeasurementUnit(r.getManipulator, body, xyz, rpy);
      
      % be able to draw stuff!
      checkDependency('lcmgl');
      obj.lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), '¯\_(ツ)_/¯');
    end
    
    function [integration_constraint, integration_equality, integration_sense] = generate_integration(obj, n, nX)
      nX_per = 12;
      integration_sense = '=';
      if (n~=1)
        % q_{k-1} - q_{k} + h\dot{q}_{k} = 0
        % velocity integrates into position for each time point
        integration_constraint = zeros(6, nX);
        integration_equality = zeros(6, 1);
        offset = nX_per*(n-2);
        integration_constraint(1:6, 1+offset:6+offset) = eye(6); %q_n
        integration_constraint(1:6, nX_per+1+offset:nX_per+6+offset) = -eye(6); %-q_n+1
        integration_constraint(1:6, nX_per+7+offset:nX_per+12+offset) = obj.dt*eye(6);%+h*qd_n+1 (backwards integration here)
      else
        warning('this shouldnt be called like this?');
        integration_constraint = zeros(0, 6);
        integration_equality = zeros(0, 6);
      end
    end
    
    function [manip_constraint, manip_equality] = generate_manip(obj, n, nX)
      state_hist = obj.state_history.getData;
      nX_per = 12;
      if (n~=1)
        manip_constraint = zeros(6, nX);
        manip_equality = zeros(6, 1);
        [Hnext,Cnext,Bnext] = obj.r.manipulatorDynamics(state_hist(1:6, n), state_hist(7:end, n));
        offset = nX_per*(n-2);
        % Hnext * (qd_n - qd_{n-1}) +
        % dt*(C_next+Gnext-B_next*u_k-J_k.'*contact_forces)
        % So dropping some terms for now (no u, no contact forces)
        manip_constraint(:, nX_per+7+offset:nX_per+12+offset) = Hnext;
        manip_constraint(:, 7+offset:12+offset) = -Hnext;
        % -J.' * contact forces
        %manip_constraint(:, nX_per+13+offset:nX_per+12+m+offset) = -J.';
        manip_equality = manip_equality - obj.dt * Cnext; % - obj.dt*Gnext;
      else
        warning('this shouldnt be called like this?');
        manip_constraint = zeros(0, nX);
        manip_equality = zeros(0, 1);
      end
    end
    
    function [imu_constraint, imu_equality] = generate_imu(obj, n, nX, x0, x1)
      nX_per = 12;
      imu_hist = obj.imu_history.getData;
      if (n~=1)
        z_k = imu_hist(:, n);
        imu_constraint = zeros(6, nX);
        %imu_equality = zeros(6, 1);
        % Taylor-expanding the measurement eq h(x) at the pair x0, x1:
        % h(xk-1, xk) = z_{k} = J(x0,x1)*([xk-1;xk]-[x0;x1]) + h(x0, x1)
        % z_{k} - h(x0,x1) + J(x0,x1)*[x0;x1] = J(x0,x1)*[xk-1;xk]
        
        % So calculate the imu jacobian numerically
        x_ext = [x0; x1];
        J = zeros(6, length(x_ext));
        out_norm = obj.imu.output(obj.r,0,0,[x1; x1; x0],[]);
        eps = 1E-6;
        add = x_ext*0;
        for j=1:length(x_ext)
          add = add*0;
          add(j) = eps;
          new_vec = x_ext + add;
          out_this = (obj.imu.output(obj.r,0,0,[new_vec(13:end); new_vec(13:end); new_vec(1:12)], []) - out_norm)/eps;
          J(:, j) = out_this;
        end
        
        % And write in the above constraint
        offset = nX_per*(n-2);
        imu_constraint(:, 1+offset:24+offset) = J;
        imu_equality = z_k - out_norm + J*[x0; x1];
      else
        warning('this shouldnt be called like this?');
        imu_constraint = zeros(0, nX);
        imu_equality = zeros(0, 1);
      end
    end
    
    function [nonpenetration_constraint, nonpenetration_equality, ...
        nonpenetration_sense] = generate_nonpenetration(obj, n, nX, m, x0)
      nX_per = 12;
      nonpenetration_constraint = zeros(m, nX);
      %nonpenetration_equality = zeros(m, 1);
      % Taylor-expanding the gap functino p(x) at the point x0:
      % p(xk) >= 0 -> J(x0)*([xk]-[x0]) + p(x0) >= 0

      % So calculate the imu jacobian numerically
      J = zeros(m, length(x0));
      q0 = x0(1:6);
      out_norm = obj.r.contactConstraints(q0,false);
      eps = 1E-6;
      add = q0*0;
      for j=1:length(q0)
        add = add*0;
        add(j) = eps;
        out_this = (obj.r.contactConstraints(q0 + add,false) - out_norm)/eps;
        J(:, j) = out_this;
      end
      % And write in the above constraint
      offset = nX_per*(n-1);
      nonpenetration_constraint(:, 1+offset:12+offset) = J;
      nonpenetration_equality = -out_norm + J*x0;
      nonpenetration_sense = '>';
    end
    
    function do_state_est(obj, t, regenerate)
      state_hist = obj.state_history.getData();
      best_guess = obj.last_sol.getData();
      %if (isempty(best_guess))
      %  best_guess = obj.init_cond;
      %else
      %  best_guess(:, 1) = obj.init_cond;
      %end
      best_guess = padarray(best_guess, size(state_hist)-size(best_guess),'replicate', 'post');
      
      N = size(state_hist, 2);
      m = obj.r.getNumContactPairs();
      
      % # of decision vars
      % x, xdot per timestep
      nX_per = 12;
      nX = N*nX_per;
      
      % # of cost terms
      % 6 manip eqs, 6 measurement eqs per pair of points
      % (measurement is a time-differencing IMU, so I'm setting this up
      % by pairs instead of per point.
      nManipCost = (N-1)*6;
      F_manip = sparse(nManipCost, nX);
      g_manip = zeros(nManipCost, 1);
      nMeasCost = (N-1)*6;
      F_meas = sparse(nMeasCost, nX);
      g_meas = zeros(nMeasCost, 1);

      % # of constraints
      % 6 integration constraints for each pair of time steps
      % plus 12 initial condition constraints
      nCons = (N-1)*6;
      A_int = sparse(nCons, nX);
      b_int = zeros(nCons, 1);
      sense_int = zeros(nCons, 1);
      
      % plus m non-penetration constraints for ALL time steps
      % (note: eventually shouldn't this be all pairwise contacts,
      % so m^2-m of these or something like that? There will be 
      % discontinuity in the gap function in general if the
      % closest-to-touching pair of bodies changes)
      nNonPenetration = N*m;
      A_nonpenetration = sparse(nNonPenetration, nX);
      b_nonpenetration = zeros(nNonPenetration, 1);
      sense_nonpenetration = zeros(nNonPenetration, 1);

      % initial condition
      A_init = padarray(obj.A_init, [12 nX]-size(obj.A_init), 0, 'post');
      b_init = obj.b_init;
      
      if (~regenerate)
        startN = N;
        F_manip(1:end-6, 1:end-nX_per) = obj.F_manip_cache.getData();
        g_manip(1:end-6) = obj.g_manip_cache.getData();
        F_meas(1:end-6, 1:end-nX_per) = obj.F_meas_cache.getData();
        g_meas(1:end-6) = obj.g_meas_cache.getData();
        A_int(1:end-6, 1:end-nX_per) = obj.A_int_cache.getData();
        b_int(1:end-6) = obj.b_int_cache.getData();
        sense_int(1:end-6) = obj.sense_int_cache.getData();
        A_nonpenetration(1:end-m, 1:end-nX_per) = obj.A_nonpenetration_cache.getData();
        b_nonpenetration(1:end-m) = obj.b_nonpenetration_cache.getData();
        sense_nonpenetration(1:end-m) = obj.sense_nonpenetration_cache.getData();
      else 
        startN = 1;
      end
      for n=startN:N
        if (nNonPenetration > 0)
          [nonpenetration_constraint, nonpenetration_equality, nonpenetration_sense] = ...
            generate_nonpenetration(obj, n, nX, m, best_guess(:, n));
          A_nonpenetration(1 + (n-1)*m : (n)*m, :) = nonpenetration_constraint;
          b_nonpenetration(1 + (n-1)*m : (n)*m) = nonpenetration_equality;
          sense_nonpenetration(1 + (n-1)*m : (n)*m) = nonpenetration_sense;
        end
        if (n ~= 1)
          % constraints
          [integration_constraint, integration_equality, integration_sense] = ...
            generate_integration(obj, n, nX);
          A_int(1 + (n-2)*6 : (n-1)*6, :) = integration_constraint;
          b_int(1 + (n-2)*6 : (n-1)*6) = integration_equality;
          sense_int(1 + (n-2)*6 : (n-1)*6) = integration_sense;
          
          % costs
          [manip_constraint, manip_equality] = ...
            generate_manip(obj, n, nX);
          F_manip(1 + (n-2)*6 : (n-1)*6, :) = manip_constraint;
          g_manip(1 + (n-2)*6 : (n-1)*6) = manip_equality;
          
          [imu_constraint, imu_equality] = ...
            generate_imu(obj, n, nX, best_guess(:, n-1), best_guess(:, n));
          F_meas(1 + (n-2)*6 : (n-1)*6, :) = imu_constraint;
          g_meas(1 + (n-2)*6 : (n-1)*6) = imu_equality;
        end
      end
      
      obj.F_manip_cache.setData(F_manip);
      obj.g_manip_cache.setData(g_manip);
      obj.F_meas_cache.setData(F_meas);
      obj.g_meas_cache.setData(g_meas);
      obj.A_int_cache.setData(A_int);
      obj.b_int_cache.setData(b_int);
      obj.sense_int_cache.setData(sense_int);
      obj.A_nonpenetration_cache.setData(A_nonpenetration);
      obj.b_nonpenetration_cache.setData(b_nonpenetration);
      obj.sense_nonpenetration_cache.setData(sense_nonpenetration);
      
      % Solve!
      % Constraint
      model.A = [A_init; A_int; A_nonpenetration];
      model.rhs = [b_init; b_int; b_nonpenetration];
      model.sense = char([obj.sense_init; sense_int; sense_nonpenetration]);
      
      % Cost -- turn equality constraint into quadratic term
      % by Q = 0.5*F.'F, c = -F.'*g
      % see http://math.stackexchange.com/questions/869204/
      % are-constrained-linear-least-squares-and-quadratic-programming-the-same-thin
      F = [F_manip; F_meas*0.00]; 
      g = [g_manip; g_meas*0.00];
      model.Q = 0.5*(F.'*F);
      model.obj = -F.'*g;
      model.objcon = g.'*g;
      
      % All terms unbounded
      model.lb = -Inf(nX, 1);
      model.ub = Inf(nX, 1);
      
      best_guess_nan_current_step = best_guess;
      best_guess_nan_current_step(:, end) = NaN;
      model.start = reshape(best_guess_nan_current_step, [numel(best_guess) 1]);
      params.outputflag = 0;
      result = gurobi(model, params)
      
      state_guess = result.x;
      
      obj.last_sol.setData(reshape(state_guess, 12, N));
      obj.full_state_est_tt.setData([obj.full_state_est_tt.getData() t]);
      % Drawing
      if ~isempty(obj.lcmgl)
        obj.lcmgl.glColor3f(1, 0, 0);
        obj.lcmgl.glPushMatrix();
        %obj.lcmgl.glLoadIdentity();
        obj.lcmgl.glTranslated(state_guess(end-11), state_guess(end-10), state_guess(end-9));
        obj.lcmgl.glRotated(state_guess(end-6)*180/pi, 0, 0, 1);
        obj.lcmgl.glRotated(state_guess(end-7)*180/pi, 0, 1, 0);
        obj.lcmgl.glRotated(state_guess(end-8)*180/pi, 1, 0, 0);
        
        obj.lcmgl.box([0,0,0], [2, 0.5, 1.0]);
        obj.lcmgl.glPopMatrix();
        %          obj.lcmgl.line3(origin(1), origin(2), origin(3), point(1,i), point(2,i), point(3,i));
        obj.lcmgl.switchBuffers;
        
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
        regenerate = mod(length(obj.t_history.getData), 10)==0 || t<0.2;
        obj.do_state_est( t, regenerate );
      end
      last_sol = obj.last_sol.getData();
      if (isempty(last_sol))
        last_sol = zeros(length(obj.getOutputFrame.getFrameByNum(3).getCoordinateNames), 1);
      else
        last_sol = full(last_sol(end-11:end));
      end
      varargout = [varargin last_sol];
      
      % benchmarking
      tim = toc;
      fprintf('Iter at %f took %f and claims yaw %f\n', t, tim, last_sol(6));
      obj.eval_time.setData([obj.eval_time.getData; tim]);
    end
  end
  
end
