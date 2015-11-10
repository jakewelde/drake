classdef ContactImplicitLogForceFixedPointProgram < NonlinearProgram

  properties
    plant % the RBM
    
    options % options, yup
    q_inds  % n x N indices for position-only state (no velocity)
    u_inds  % m x N indices for time
    dynamic_constraints = {};
    constraints = {};
    
    nC
    nD % number of friction elements per contact
    
    l_inds % orderered [lambda_N;lambda_f1;lambda_f2] for each contact sequentially
           % in this variant, represent LOG of contact forces
    lfi_inds % nD x nC indexes into lambda for each time step
    lambda_mult
    ljl_inds  % joint limit forces
    jl_lb_ind  % joint indices where the lower bound is finite
    jl_ub_ind % joint indices where the lower bound is finite
    nJL % number of joint limits = length([jl_lb_ind;jl_ub_ind])
    
    nonlincompl_constraint
    nonlincompl_slack_inds
    
    jlcompl_constraint
    jlcompl_slack_inds
    
    contact_epsilon = 100;
    
  end
  
  methods
    function obj = ContactImplicitLogForceFixedPointProgram(plant,x_dimensions_to_ignore, options)
      % @param sys an rbm
      % @param x_dimensions_to_ignore if this is specified then xdot need
      % not be zero in the ignored dimensions (e.g. useful for finding a
      % trim condition of an aircraft instead of a true fixed point)
      if nargin<3, options=struct(); end
      
      if ~isfield(options,'nlcc_mode')
        options.nlcc_mode = 3;
      end
      if ~isfield(options,'lincc_mode')
        options.lincc_mode = 1;
      end
      if ~isfield(options,'compl_slack')
        options.compl_slack = 0;
      end
      if ~isfield(options,'lincompl_slack')
        options.lincompl_slack = 0;
      end
      if ~isfield(options,'jlcompl_slack')
        options.jlcompl_slack = 0;
      end
      if ~isfield(options,'lambda_mult')
        options.lambda_mult = 1;
      end
      if ~isfield(options,'lambda_jl_mult')
        options.lambda_jl_mult = 1;
      end
      if ~isfield(options,'active_collision_options')
        options.active_collision_options.terrain_only = false;
      end
      if ~isfield(options, 'multiple_contacts')
        options.multiple_contacts = false;
      end

      typecheck(plant,'RigidBodyManipulator');
      if ~isTI(plant), error('only makes sense for time invariant systems'); end
      if nargin<2
        x_dimensions_to_ignore = [];
      else
        'this does nothing'
      end
      
      nQ = getNumPositions(plant);
      nU = getNumInputs(plant);

      positionNames = getCoordinateNames(plant.getStateFrame);
      positionNames = positionNames(1:nQ);
      obj = obj@NonlinearProgram(nQ+nU,vertcat(positionNames, getCoordinateNames(plant.getInputFrame)));

      obj.options = options;
      obj.plant = plant;

      obj.q_inds = reshape((1:nQ),nQ,1);
      obj.u_inds = reshape(nQ + (1:nU),nU,1);

      obj.nC = obj.plant.getNumContactPairs;
      [~,normal,d] = obj.plant.contactConstraints(zeros(obj.plant.getNumPositions,1), obj.options.multiple_contacts);
      obj.nD = 2*length(d);
      assert(size(normal,2) == obj.nC); % just a double check
      
      % contact forces along friction cone bases
      nContactForces = obj.nC*(1 + obj.nD);    
      obj.l_inds = reshape(obj.num_vars + (1:nContactForces),nContactForces,1);
      obj = obj.addDecisionVariable(nContactForces);
      
      %???
      obj.lfi_inds = zeros(obj.nD,obj.nC);
      for i=1:obj.nC,
        obj.lfi_inds(:,i) = (2:1+obj.nD)' + (i-1)*(1+obj.nD)*ones(obj.nD,1);
      end
      
      obj.nJL = obj.plant.getNumJointLimitConstraints();
      obj.ljl_inds = reshape(obj.num_vars + (1:obj.nJL),obj.nJL,1);
      
      % joint limit constraints
      [jl_lb,jl_ub] = obj.plant.getJointLimits();
      obj.jl_lb_ind = find(jl_lb ~= -inf);
      obj.jl_ub_ind = find(jl_ub ~= inf);
      obj = obj.addDecisionVariable(obj.nJL);

      obj = addDynamicConstraints(obj,x_dimensions_to_ignore);
      
      % todo: state and input constraints
    end
    
    function obj = addDynamicConstraints(obj,x_dimensions_to_ignore)
      nQ = obj.plant.getNumPositions();
      nU = obj.plant.getNumInputs();

      %      state input   contact forces      joint limits
      n_vars = nQ + nU + obj.nC*(1+obj.nD) + obj.nJL;
      dynamicConstraint = FunctionHandleConstraint(0*ones(nQ,1),0*ones(nQ,1),n_vars,@obj.dynamics_constraint_fun);
      dynInds = {obj.q_inds;obj.u_inds;obj.l_inds;obj.ljl_inds};
      obj = obj.addConstraint(dynamicConstraint, dynInds);
      
      [~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(zeros(nQ,1),obj.options.multiple_contacts,obj.options.active_collision_options);
      
      if obj.nC > 0
        lambda_N_inds = obj.l_inds(1:(obj.nD+1):end);
        
        nonpen = FunctionHandleConstraint(0*ones(obj.nC, 1), Inf*ones(obj.nC, 1), nQ, @nonpen_fun);
        obj = obj.addConstraint(nonpen,[obj.q_inds]);
        
        nonpen_force_compl = FunctionHandleConstraint(-Inf*ones(obj.nC, 1), log(obj.contact_epsilon)*ones(obj.nC, 1), nQ+obj.nC, @nonpen_force_compl_fun);
        obj = obj.addConstraint(nonpen_force_compl,[obj.q_inds;lambda_N_inds]);
        
        % lambda_fi >= 0
        %friction_cone_constraint = FunctionHandleConstraint(0*ones(obj.nC,1), Inf*ones(obj.nC,1), obj.nC*(obj.nD+1), @friction_cone_fun);
        %obj = obj.addConstraint(friction_cone_constraint,[obj.l_inds]);
        
        obj = obj.addConstraint(BoundingBoxConstraint(-Inf*ones(obj.nC*(obj.nD),1), log(1E-3)*ones(obj.nC*(obj.nD),1)), [obj.lfi_inds]);
      end
    
      if obj.nJL > 0
        % joint limit linear complementarity constraint
        % lambda_jl /perp [q - lb_jl; -q + ub_jl]
        W_jl = zeros(obj.nJL);
        [r_jl,M_jl] = jointLimitConstraints(obj.plant,zeros(nQ,1));
        obj.jlcompl_constraint = LinearComplementarityConstraint(W_jl,r_jl,M_jl,obj.options.lincc_mode,obj.options.jlcompl_slack);
        obj.jlcompl_slack_inds = obj.num_vars+1:obj.num_vars + obj.jlcompl_constraint.n_slack;
        obj = obj.addConstraint(obj.jlcompl_constraint,[obj.q_inds;obj.ljl_inds]);
      end
     
      % nonpenetration:
      %   phi(q) >= 0
      % x = [q]
      function [f,df] = nonpen_fun(q)
        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,obj.options.multiple_contacts,obj.options.active_collision_options);

        f = phi;
        df = n;
      end
      
      % nonlinear complementarity constraints:
      %   e^lambda_N .* phi(q) = 0
      %   --> lambda_N + ln(phi(q)) <= ln(eps)
      % x = [q, lambda_N]
      function [f,df] = nonpen_force_compl_fun(y)
        nq = obj.plant.getNumPositions;
        q = y(1:nq);
        lambda_n = y(nq+1:end);

        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,obj.options.multiple_contacts,obj.options.active_collision_options);

        f = zeros(obj.nC,1);
        df = zeros(obj.nC,nq+obj.nC);

        f(:,1) = lambda_n + log(phi);
        df(:,1:nq) = repmat(1./phi, 1, nq) .* n;
        df(:,nq+1:nq+obj.nC) = eye(obj.nC);
      end
      
      % friction cone constraints
      % mu*lambda_N >= sum(lambda_fi)
      % --> mu*exp(lambda_N) >= sum(exp(lambda_fi))
      function [f,df] = friction_cone_fun(y)
        nq = obj.plant.getNumPositions;
        x = y(1:nq+obj.nC);
        
        f = zeros(obj.nC, 1);
        df = zeros(obj.nC, obj.nC*(obj.nD+1));
        
        for i=1:obj.nC
          offset = (i-1)*(obj.nD+1);
           f(i) = mu(i) * exp(y(offset + 1)) - sum(exp(y(offset + [2:(obj.nD+1)])));
           df(i, offset+1:offset+obj.nD+1) = [mu(i); -ones(obj.nD, 1)];
        end
      end
      
    end

    function [f,df] = dynamics_constraint_fun(obj,q0,u,lambda,lambda_jl)
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nu = obj.plant.getNumInputs;
      nl = length(lambda);
      njl = length(lambda_jl);
      
      lambda = exp(lambda)*obj.options.lambda_mult;
      lambda_jl = lambda_jl*obj.options.lambda_jl_mult;
        
      assert(nq == nv) % not quite ready for the alternative
      
      v0 = q0*0;
      
      [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(q0, v0);
      
      BuminusC = B*u-C; % functionally B*u-G, since v0 = 0

      if nu>0,  
        dBuminusC = matGradMult(dB,u) - dC;
      else
        dBuminusC = -dC;
      end
      
      % 0 = (B*u - C) + n^T lambda_N + d^T * lambda_f
      fv = BuminusC;
      % [q0 v0 u l ljl]
      
      dfv = [zeros(nq,nq), B, zeros(nv,nl+njl)] + ...
        [dBuminusC(:, 1:nq) zeros(nq,nu+nl+njl)];
      
      if nl>0
        [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q0,obj.options.multiple_contacts,obj.options.active_collision_options);
        % construct J and dJ from n,D,dn, and dD so they relate to the
        % lambda vector
        J = zeros(nl,nq);
        J(1:1+obj.nD:end,:) = n;
        dJ = zeros(nl*nq,nq);
        dJ(1:1+obj.nD:end,:) = dn;
        
        for j=1:length(D),
          J(1+j:1+obj.nD:end,:) = D{j};
          dJ(1+j:1+obj.nD:end,:) = dD{j};
        end

        fv = fv + J'*lambda;
        dfv(:,1:nq) = dfv(:,1:nq) + matGradMult(dJ,lambda,true);
        dfv(:,1+nq+nu:nq+nu+nl) = J'*obj.options.lambda_mult;
      end
      
      if njl>0
        [~,J_jl,dJ_jl] = jointLimitConstraints(obj.plant,q0);
        
        fv = fv + J_jl'*lambda_jl;
        dfv(:,1+nq+nu+nl:nq+nu+nl+njl) = J_jl'*obj.options.lambda_jl_mult;
        %necessary? don't think so, since dJ_jl is zero unless we're at
        % the constraint boundary, at which point the complementarity takes
        % over
        %dfv(:,1:nq) = dfv(:,1:nq) + reshape(dJ_jl'*lambda_jl, 36, 36);
      end
      
      % Update dfv for lambda to account for the exponentiation we did
      % earlier. df/dlambda = df/dexp(lambda) dexp(lambda)/dlambda
      % = dfv * e^(lambda)
      % (lambda has been self-assigned e^lambda so just mult by it again)
      dfv(:, 1+nq+nu:nq+nu+nl) = dfv(:, 1+nq+nu:nq+nu+nl) * repmat(lambda, 1, nl);
      f = [fv];
      df = [dfv];
    end
    
    
    % 
    function z0 = getInitialVars(obj,stateguess)
      z0 = zeros(obj.num_vars,1);
      
      if ~isfield(stateguess,'q')
        nQ = getNumPositions(obj.plant);
        stateguess.q = 0.01*randn(nQ,1);
      end
      z0(obj.q_inds) = stateguess.q;
      
      kinsol = obj.plant.doKinematics(stateguess.q);
      [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(stateguess.q,obj.options.multiple_contacts,obj.options.active_collision_options);
      [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(stateguess.q, 0*stateguess.q);
      nU = getNumInputs(obj.plant);
      nl = numel(obj.l_inds);
      nq = obj.plant.getNumPositions;
      if (numel(phi) > 0)
        J = zeros(nl,nq);
        J(1:1+obj.nD:end,:) = n;
        for j=1:length(D),
          J(1+j:1+obj.nD:end,:) = D{j};
        end
      end
      
      if ~isfield(stateguess,'l') && ~isfield(stateguess, 'u')
          % jointly calculate contact force and control input guess
          lu = [J.' B] \ C;
          stateguess.l = log(lu(1:nl)+1E-12);
          stateguess.u = lu((nl+1):end);
      end
      
      if obj.nC > 0
        if ~isfield(stateguess,'l')
          % calculate an l guess: zero where phi!=0,
          % solution to manip eq otherwise
          %stateguess.l = log(1E-9)*ones(numel(obj.l_inds),1);
          stateguess.l = log((J.') \ (C - B*zeros(nU,1)));
        end
        z0(obj.l_inds) = stateguess.l;
      end
      if obj.nJL > 0
        if ~isfield(stateguess,'ljl')
          % guess that we start out away from joint lims
          stateguess.ljl = 0;
        end
        z0(obj.ljl_inds) = stateguess.ljl;
      end
      
      % explain as much of the rest as possible with joint torques
      if ~isfield(stateguess,'u')
        if (nU > 0)
          stateguess.u = B \ (C - J.'*exp(stateguess.l));
        else
          stateguess.u = 0;
        end
      end
      z0(obj.u_inds) = stateguess.u;    
    end

    
    function [qstar,ustar,lstar,info,F] = findFixedPoint(obj,guess,v)
      if nargin<3, v = []; end
      if nargin<2, guess = struct(); end

      if ~isempty(v)
        prog = obj.addDisplayFunction(@(x)v.drawWrapper(0,[x;zeros(obj.plant.getNumPositions, 1)]));
      end

      z0 = prog.getInitialVars(guess);
      [z,F,info,infeasible_constraint_name] = prog.solve(z0);
      qstar = z(obj.q_inds);
      ustar = Point(obj.plant.getInputFrame,z(obj.u_inds));
      lstar = z(obj.l_inds);
    end

    function obj = addStateCost(obj,state_cost_function)
      nQ = obj.plant.getNumPositions();
      state_cost = FunctionHandleObjective(nQ,state_cost_function, -1);
      %state_cost.grad_method = 'numerical';
      obj = obj.addCost(state_cost,{obj.q_inds});
    end

    function obj = addInputCost(obj,input_cost_function)
      %state_cost.grad_method = 'numerical';
      obj = obj.addCost(input_cost_function,{obj.u_inds});
    end
    
  end

end