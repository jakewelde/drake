classdef ContactImplicitFixedPointUnconstrainedProgram < NonlinearProgram

  properties
    plant % the RBM
    
    options % options, yup
    q_inds  % n x N indices for position-only state (no velocity)
    u_inds  % m x N indices for time
    dynamic_constraints = {};
    constraints = {};
    
    nC
    nD % number of friction elements per contact
    
    l_inds % orderered [lambda_N;lambda_f1;lambda_f2;...] for each contact sequentially
    lfi_inds % nD x nC indexes into lambda for each time step
    lambda_mult
    ljl_inds  % joint limit forces
    jl_lb_ind  % joint indices where the lower bound is finite
    jl_ub_ind % joint indices where the lower bound is finite
    nJL % number of joint limits = length([jl_lb_ind;jl_ub_ind])
    
    total_scale = 1;
    dynamics_scale = 1.0;
    nonlincompl_scale = 100.0;
    lincompl_scale = 100.0;
    penetration_scale = 10.0;
    jlcompl_scale = 100.0;
    
    
  end
  
  methods
    function obj = ContactImplicitFixedPointUnconstrainedProgram(plant,x_dimensions_to_ignore, options)
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
      if ~isfield(options, 'scaling')
          options.scaling = 1;
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

      scale = obj.total_scale * options.scaling;
      obj.dynamics_scale = obj.dynamics_scale*scale;
      obj.nonlincompl_scale = obj.nonlincompl_scale*scale;
      obj.lincompl_scale = obj.lincompl_scale*scale;
      obj.penetration_scale = obj.penetration_scale*scale;
      
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
      
      % separate out Lambda_fi inds
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
      
      dynInds = {obj.q_inds;obj.u_inds;obj.l_inds;obj.ljl_inds};
      dynamicCost = FunctionHandleConstraint(0,0,n_vars,@obj.dynamics_constraint_fun);
      %dynamicCost.grad_method = 'numerical';
      obj = obj.addCost(dynamicCost, dynInds);
      
      [~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(zeros(nQ,1),obj.options.multiple_contacts,obj.options.active_collision_options);
      
      if obj.nC > 0
        lambda_inds = obj.l_inds;
          
        % positivity of contact forces
        bounds = ones(numel(obj.l_inds),1);
        obj = obj.addConstraint(BoundingBoxConstraint(bounds*0, bounds*Inf), obj.l_inds);

        nonlincomplCost = FunctionHandleConstraint(0,0,nQ + obj.nC*(1+obj.nD),@obj.nonlincompl_fun);
        %nonlincomplCost.grad_method = 'numerical';
        obj = obj.addCost(nonlincomplCost,{obj.q_inds;obj.l_inds});
        
        % friction cone description
        %  0 <= mu*lambda_N - sum(lambda_fi)
        A = zeros(obj.nC, obj.nC*(obj.nD+1));
        for k=1:obj.nC
           A(k, (1:(obj.nD+1)) + ones(1,obj.nD+1)*(k-1)*(obj.nD+1)) = ...
               [mu(k), -ones(1, obj.nD)];
        end
        obj = obj.addConstraint(LinearConstraint(zeros(obj.nC,1), ones(obj.nC,1)*Inf, A), obj.l_inds);
      end
    
      if obj.nJL > 0
        % joint limit linear complementarity constraint
        % lambda_jl /perp [q - lb_jl; -q + ub_jl]
        % positivity of joint limits forces
        obj = obj.addConstraint(BoundingBoxConstraint(zeros(obj.nJL,1), ones(obj.nJL,1)*Inf), [obj.ljl_inds]);
        % joint limits themselves
        [r_jl, M_jl] = jointLimitConstraints(obj.plant,zeros(nQ,1));
        obj = obj.addConstraint(LinearConstraint(-r_jl(1:obj.nJL/2), r_jl((obj.nJL/2+1):end), M_jl(1:obj.nJL/2,:)), [obj.q_inds]);
        % restoring force at joint limits with relaxed complementarity
        jlcomplCost = FunctionHandleConstraint(0,0,nQ + obj.nJL, @obj.jlcompl_fun);
        obj = obj.addCost(jlcomplCost, {obj.q_inds; obj.ljl_inds});
      end
      
    end
    
    % joint limit complementarity constraints:
    %   lambda_jl /perp [q - lb_jl; -q + ub_jl]
    function [f,df] = jlcompl_fun(obj, q, lambda)
        nq = obj.plant.getNumPositions;
        v = q*0;
        
        [r_jl,M_jl] = jointLimitConstraints(obj.plant,zeros(nq,1));
        
        f = 0;
        df = zeros(1,nq+obj.nJL);
        
        err = lambda .* (M_jl*q + r_jl);
        
        f = f + err.'*err;
        df(1:nq) = df(1:nq) + 2 * err.' * (repmat(lambda, 1, nq) .* M_jl);
        df(nq+1:nq+obj.nJL) = 2 * err.' * repmat(M_jl*q + r_jl, 1, obj.nJL);
        
        f = f * obj.jlcompl_scale;
        df = df * obj.jlcompl_scale;
    end
    
    % nonlinear complementarity constraints:
    %   lambda_N /perp phi(q) -> lambda_N * phi(q) = 0, penalized quadratically
    % x = [q;lambda_N;lambda_F1;lambda_f2](each contact sequentially)
    function [f,df] = nonlincompl_fun(obj, q, lambda)
        nq = obj.plant.getNumPositions;
        v = q*0;
        
        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,obj.options.multiple_contacts,obj.options.active_collision_options);
        
        f = 0;
        df = zeros(1,nq+obj.nC*(1+obj.nD));
        
        % nonpenetration
        err = phi(phi<0);
        f = f + obj.penetration_scale*(err.'*err);
        df(1:nq) = df(1:nq) + sum(2*obj.penetration_scale*repmat(err,1,nq).*n(phi<0,:));
        
        % lambda_N * phi(q) = 0
        err = phi.*lambda(1:obj.nD+1:end);
        
        % L2
        f = f + err.'*err * obj.nonlincompl_scale;
        df(1:nq) = df(1:nq) + obj.nonlincompl_scale*sum( 2*repmat(err, 1, nq).*n.*repmat(lambda(1:obj.nD+1:end), 1, nq) );
        lambdaNInds = (nq+1):obj.nD+1:(nq+obj.nC*(obj.nD+1));
        df(lambdaNInds) = df(lambdaNInds) + ...
            obj.nonlincompl_scale*(2*err.*phi).';
        
        f = f;
        df = df;
    end
    
    function [f,df] = dynamics_constraint_fun(obj,q0,u,lambda,lambda_jl)
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nu = obj.plant.getNumInputs;
      nl = length(lambda);
      njl = length(lambda_jl);
      
      lambda = lambda*obj.options.lambda_mult;
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
        [~,J_jl] = jointLimitConstraints(obj.plant,q0);
        
        fv = fv + J_jl'*lambda_jl;
        dfv(:,1+nq+nu+nl:nq+nu+nl+njl) = J_jl'*obj.options.lambda_jl_mult;
      end
      
      % apply norm
      f = fv.'*fv;
      df = (2*dfv.'*fv).';
      
      f = f * obj.dynamics_scale;
      df = df * obj.dynamics_scale;
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
          % todo: joint limits might be needed here.
          lu = zeros(nU+nl, 1);
          legal_J = J(phi==0, :);
          lu([find(phi==0); [(nl+1):nl+nU]]) = [legal_J.' B] \ C;
          stateguess.l = lu(1:nl);
          stateguess.u = lu((nl+1):end);
      end
      
      if obj.nC > 0
        if ~isfield(stateguess,'l')
          % calculate an l guess: zero where phi!=0,
          % solution to manip eq otherwise
          stateguess.l = zeros(numel(obj.l_inds),1);
          stateguess.l = (J.') \ (C - B*zeros(nU,1));
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
          stateguess.u = B \ (C - J.'*stateguess.l);
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
%       %z0 = z0 + rand;
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