classdef ContactImplicitFixedPointProgram < NonlinearProgram

  properties
    plant % the RBM
    
    options % options, yup
    q_inds  % n x N indices for position-only state (no velocity)
    u_inds  % m x N indices for time
    dynamic_constraints = {};
    constraints = {};
    
    nC
    nD % number of friction elements per contact
    
    l_inds % orderered [lambda_N;lambda_f1;lambda_f2;...;gamma] for each contact sequentially
    lfi_inds % nD x nC indexes into lambda for each time step
    lambda_mult
    ljl_inds  % joint limit forces
    jl_lb_ind  % joint indices where the lower bound is finite
    jl_ub_ind % joint indices where the lower bound is finite
    nJL % number of joint limits = length([jl_lb_ind;jl_ub_ind])
    
    nonlincompl_constraint
    nonlincompl_slack_inds
    
  end
  
  methods
    function obj = ContactImplicitFixedPointProgram(plant,x_dimensions_to_ignore, options)
      % @param sys an rbm
      % @param x_dimensions_to_ignore if this is specified then xdot need
      % not be zero in the ignored dimensions (e.g. useful for finding a
      % trim condition of an aircraft instead of a true fixed point)
      if nargin<3, options=struct(); end
      
      if ~isfield(options,'nlcc_mode')
        options.nlcc_mode = 2;
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
      nContactForces = obj.nC*(2 + obj.nD);    
      obj.l_inds = reshape(obj.num_vars + (1:nContactForces),nContactForces,1);
      obj = obj.addDecisionVariable(nContactForces);
      
      %???
      obj.lfi_inds = zeros(obj.nD,obj.nC);
      for i=1:obj.nC,
        obj.lfi_inds(:,i) = (2:1+obj.nD)' + (i-1)*(2+obj.nD)*ones(obj.nD,1);
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

      nq = obj.plant.getNumPositions();
      fixedPointConstraint = ConstantConstraint(zeros(nq, 1));
      obj.addConstraint(fixedPointConstraint, obj.q_inds);
    end
    
    function obj = addDynamicConstraints(obj,x_dimensions_to_ignore)
      nQ = obj.plant.getNumPositions();
      nU = obj.plant.getNumInputs();

      %      state input   contact forces      joint limits
      n_vars = nQ + nU + obj.nC*(2+obj.nD) + obj.nJL;
      dynamicConstraint = FunctionHandleConstraint(zeros(nQ,1),zeros(nQ,1),n_vars,@obj.dynamics_constraint_fun);
      dynInds = {obj.q_inds;obj.u_inds;obj.l_inds;obj.ljl_inds};
      obj = obj.addConstraint(dynamicConstraint, dynInds);

      stabilityConstraint = FunctionHandleConstraint(-Inf, ones(1,1),n_vars,@obj.dynamics_linstability_constraint_fun);
      stabilityConstraint.grad_method = 'numerical';
      dynInds = {obj.q_inds;obj.u_inds;obj.l_inds;obj.ljl_inds};
      obj = obj.addConstraint(stabilityConstraint, dynInds);
      
      [~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(zeros(nQ,1),obj.options.multiple_contacts,obj.options.active_collision_options);
      
      if obj.nC > 0
        % indices for gamma
        gamma_inds = obj.l_inds(obj.nD+2:obj.nD+2:end,1);
        % awkward way to pull out these indices, for (i) lambda_N and
        % lambda_f
        lambda_inds = obj.l_inds(repmat((1:1+obj.nD)',obj.nC,1) + kron((0:obj.nC-1)',(2+obj.nD)*ones(obj.nD+1,1)),1);
        
        obj.nonlincompl_constraint = NonlinearComplementarityConstraint(@nonlincompl_fun,nQ + obj.nC,obj.nC*(1+obj.nD),obj.options.nlcc_mode,obj.options.compl_slack);
        obj.nonlincompl_slack_inds = obj.num_vars+1:obj.num_vars + obj.nonlincompl_constraint.n_slack;          
        obj = obj.addConstraint(obj.nonlincompl_constraint,[obj.q_inds;gamma_inds;lambda_inds]);
        
        % linear complementarity constraint
        %   gamma /perp mu*lambda_N - sum(lambda_fi)
        %
        %  Generate terms W,r,M,gamma_inds so that
        %  gamma = y(gamma_inds)
        %  Wz+Mx+r = mu*lambda_N - sum(lambda_fi)
        r = zeros(obj.nC,1);
        W = zeros(obj.nC,obj.nC);
        M = zeros(obj.nC,obj.nC*(1+obj.nD));
        for k=1:obj.nC,
          M(k,1 + (k-1)*(1+obj.nD)) = mu(k);
          M(k,(2:obj.nD+1) + (k-1)*(1+obj.nD)) = -ones(obj.nD,1);
        end
        
        lincompl_constraint = LinearComplementarityConstraint(W,r,M,obj.options.lincc_mode,obj.options.lincompl_slack);
        obj = obj.addConstraint(lincompl_constraint,[lambda_inds;gamma_inds]);
      end
    
      if obj.nJL > 0
        % joint limit linear complementarity constraint
        % lambda_jl /perp [q - lb_jl; -q + ub_jl]
        W_jl = zeros(obj.nJL);
        [r_jl,M_jl] = jointLimitConstraints(obj.plant,zeros(nq,1));
        jlcompl_constraint = LinearComplementarityConstraint(W_jl,r_jl,M_jl,obj.options.lincc_mode,obj.options.jlcompl_slack);
        
        obj = obj.addConstraint(jlcompl_constraint,[obj.q_inds;obj.ljl_inds]);
      end
      
      % nonlinear complementarity constraints:
      %   lambda_N /perp phi(q)
      %   lambda_fi /perp gamma + Di*psi(q,v)
      % x = [q;v;gamma]
      % z = [lambda_N;lambda_F1;lambda_f2] (each contact sequentially)
      function [f,df] = nonlincompl_fun(y)
        nq = obj.plant.getNumPositions;
        nv = obj.plant.getNumVelocities;
        x = y(1:nq+obj.nC);
        z = y(nq+obj.nC+1:end);
        gamma = x(nq+1:end);
        
        q = x(1:nq);
        v = q*0;
        
        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,obj.options.multiple_contacts,obj.options.active_collision_options);
        
        f = zeros(obj.nC*(1+obj.nD),1);
        df = zeros(obj.nC*(1+obj.nD),nq+obj.nC*(2+obj.nD));
        
        f(1:1+obj.nD:end) = phi;
        df(1:1+obj.nD:end,1:nq) = n;
        for j=1:obj.nD,
          f(1+j:1+obj.nD:end) = gamma+D{j}*v;
          df(1+j:1+obj.nD:end,nq+(1:obj.nC)) = eye(size(D{j},1));  %d/dgamma
          df(1+j:1+obj.nD:end,1:nq) = matGradMult(dD{j},v);%d/dq
        end
      end
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
        J(1:2+obj.nD:end,:) = n;
        dJ = zeros(nl*nq,nq);
        dJ(1:2+obj.nD:end,:) = dn;
        
        for j=1:length(D),
          J(1+j:2+obj.nD:end,:) = D{j};
          dJ(1+j:2+obj.nD:end,:) = dD{j};
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
      
      f = [fv];
      df = [dfv];
    end
    
    function f = dynamics_linstability_constraint_fun(obj,q0,u,lambda,lambda_jl)
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

      d = [];
      if nl>0
        [phi,normal,d,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q0,obj.options.multiple_contacts,obj.options.active_collision_options);
        for j=1:numel(d)
           d{j} = d{j}(:, phi<0); 
        end
        % construct J and dJ from n,D,dn, and dD so they relate to the
        % lambda vector
        J = zeros(nl,nq);
        J(1:2+obj.nD:end,:) = n;
        dJ = zeros(nl*nq,nq);
        dJ(1:2+obj.nD:end,:) = dn;

        for j=1:length(D),
          J(1+j:2+obj.nD:end,:) = D{j};
          dJ(1+j:2+obj.nD:end,:) = dD{j};
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
      
      
      % calculate qdd to help guide us in if we're not feasible yet
      %qdd = pinv(H)*fv;
      %dfv(:,1:nq) = dfv(:, 1:nq) - reshape(dH(:, 1:nq) * qdd, nq, nq);
      % qdd = A * (q - q0) linearization
      % sol like (q-q0) = c * exp(lambda*t)
      % plugging back in leads to (A - lambda^2)*c = 0
      % i.e. lambdas = sqrt'd eigenvalues of A
      % want them all nonpositive for stability
      %[V, theseeigs] = eig(pinv(H)*dfv(:, 1:nq));
      %theseeigs = diag(theseeigs);
      %theseeigs = abs(real(theseeigs.^(1/2)));
      %theseeigs = real(eig(-pinv(H)*reshape(dH(:, 1:nq), nq, nq)));
      %f = max(theseeigs)
      %f(1:numel(theseeigs)) = theseeigs;
      
      % qdd / dq
      min(phi)
      d
      Hinv = pinv(H);
      qdd_dq = -Hinv * reshape(dH(:, 1:nq) * Hinv * fv, nq, nq) + Hinv * dfv(:,1:nq);
      % qdd / du
      qdd_du = Hinv * B;
      % qdd / dlambda
      % probably not handling joint lims correctly yet
      qdd_dlambda = Hinv * J.';
      A = [qdd_dq qdd_du qdd_dlambda];
      
      [V, theseeigs] = eig(qdd_dq);
      theseeigs = diag(theseeigs);
      theseeigs = abs(real(theseeigs.^(1/2)));
      f = max(theseeigs)
      
    end
    
    % 
    function z0 = getInitialVars(obj,stateguess)
      z0 = zeros(obj.num_vars,1);

      if ~isfield(stateguess,'u')
        nU = getNumInputs(obj.plant);
        stateguess.u = 0.01*randn(nU,1);
      end
      z0(obj.u_inds) = stateguess.u;      
      
      if ~isfield(stateguess,'q')
        nQ = getNumPositions(obj.plant);
        stateguess.q = 0.01*randn(nQ,1);
      end
      z0(obj.q_inds) = stateguess.q;
      
      if obj.nC > 0
        if ~isfield(stateguess,'l')
          stateguess.l = 0;
        end
        z0(obj.l_inds) = stateguess.l;
      end
      if obj.nJL > 0
        if ~isfield(traj_init,'ljl')
          stateguess.ljl = 0;
        end
        z0(obj.ljl_inds) = stateguess.ljl;
      end
            
      if obj.nC > 0
        gamma_inds = obj.l_inds(obj.nD+2:obj.nD+2:end, 1);
        lambda_inds = obj.l_inds(repmat((1:1+obj.nD)',obj.nC,1) + kron((0:obj.nC-1)',(2+obj.nD)*ones(obj.nD+1,1)),1);          
        if ~isempty(obj.nonlincompl_slack_inds)
          z0(obj.nonlincompl_slack_inds) = obj.nonlincompl_constraint.slack_fun(z0([obj.q_inds;gamma_inds;lambda_inds]));
        end
      end
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

  end

end