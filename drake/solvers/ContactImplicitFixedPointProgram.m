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
    
    l_inds % orderered [lambda_N;lambda_f1;lambda_f2] for each contact sequentially
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
    
  end
  
  methods
    function obj = ContactImplicitFixedPointProgram(plant,x_dimensions_to_ignore, options)
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
      %dynamicCost = FunctionHandleConstraint(0,0,n_vars,@(q,u,l,ljl)norm(obj.dynamics_constraint_fun(q,u,l,ljl)));
     % dynamicCost.grad_method = 'numerical';
    %  obj = obj.addCost(dynamicCost, dynInds);
      
      
%       stabilityConstraint = FunctionHandleConstraint(-Inf, ones(1,1)*-1E-3,n_vars,@obj.dynamics_linstability_constraint_fun);
%       stabilityConstraint.grad_method = 'numerical';
%       dynInds = {obj.q_inds;obj.u_inds;obj.l_inds;obj.ljl_inds};
%       obj = obj.addConstraint(stabilityConstraint, dynInds);
      
      [~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(zeros(nQ,1),obj.options.multiple_contacts,obj.options.active_collision_options);
      
      if obj.nC > 0
        lambda_N_inds = obj.l_inds(1:(obj.nD+1):end);
        
        obj.nonlincompl_constraint = NonlinearComplementarityConstraint(@nonlincompl_fun,nQ,obj.nC,obj.options.nlcc_mode,obj.options.compl_slack);
        obj.nonlincompl_slack_inds = obj.num_vars+1:obj.num_vars + obj.nonlincompl_constraint.n_slack;
        obj = obj.addConstraint(obj.nonlincompl_constraint,[obj.q_inds;lambda_N_inds]);
        
        % lambda_fi >= 0
        obj = obj.addConstraint(BoundingBoxConstraint(zeros(obj.nC*(obj.nD),1), ones(obj.nC*(obj.nD),1)*Inf), [obj.lfi_inds]);
        
        
        % friction cone description
        %  0 <= mu*lambda_N - sum(lambda_fi)
        A = zeros(obj.nC, obj.nC*(obj.nD+1));
        for k=1:obj.nC
           A(k, (1:(obj.nD+1)) + ones(1,obj.nD+1)*(k-1)) = ...
               [mu(k), -ones(1, obj.nD)];
        end
        obj = obj.addConstraint(LinearConstraint(zeros(obj.nC,1), ones(obj.nC,1)*Inf, A), obj.l_inds);
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
      
      % nonlinear complementarity constraints:
      %   lambda_N /perp phi(q)
      % x = [q]
      % z = [lambda_N] (each contact sequentially)
      function [f,df] = nonlincompl_fun(y)
        nq = obj.plant.getNumPositions;
        nv = obj.plant.getNumVelocities;
        x = y(1:nq+obj.nC);
        z = y(nq+obj.nC+1:end);
        
        q = x(1:nq);
        v = q*0;
        
        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,obj.options.multiple_contacts,obj.options.active_collision_options);
        
        f = zeros(obj.nC,1);
        df = zeros(obj.nC,nq+obj.nC);
        
        f(:,1) = phi;
        df(:,1:nq) = n;
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
      
      f = [fv];
      df = [dfv];
    end
    
    function f = dynamics_linstability_constraint_fun(obj,q0,u,lambda,lambda_jl)
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nu = obj.plant.getNumInputs;
      nl = length(lambda);
      njl = length(lambda_jl);
      
      kinsol = obj.plant.doKinematics(q0);

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
        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q0,obj.options.multiple_contacts,obj.options.active_collision_options);
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
%       min(phi)
%       d
%       Hinv = pinv(H);
%       qdd_dq = -Hinv * reshape(dH(:, 1:nq) * Hinv * fv, nq, nq) + Hinv * dfv(:,1:nq);
%       % qdd / du
%       qdd_du = Hinv * B;
%       % qdd / dlambda
%       % probably not handling joint lims correctly yet
%       qdd_dlambda = Hinv * J.';
%       A = [qdd_dq qdd_du qdd_dlambda];
%       
%       [V, theseeigs] = eig(qdd_dq);
%       theseeigs = diag(theseeigs);
%       theseeigs = abs(real(theseeigs.^(1/2)));
%       f = max(theseeigs)

        % this doesn't work because it doesn't help guide the sol down
        % to the "good" contact configs.
        % WHYYY IS IT SO RESISTANT
        f = 0;
        for i=2:max(idxB)
           bodycom = obj.plant.forwardKin(kinsol, 1, obj.plant.body(i).com);
           good_normals_per_body = [normal(:, phi<0.001 & idxB==i) (-obj.plant.gravity/norm(obj.plant.gravity))]; 
           if (size(unique(good_normals_per_body.', 'rows'), 1) >= 4)
               K = convhulln(good_normals_per_body.');
               for j=1:size(K, 1)
                  facet = good_normals_per_body(:,K(j,:)); 
                  pnormal = cross(facet(:,1), facet(:,2));
                  if (norm(pnormal) > 0.1)
                    d = -pnormal.' * facet(:,1);
                    dist = pnormal.' * bodycom + d;
                    f = f - dist;
                  end
               end
           end
        end
      
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
        if ~isfield(stateguess,'ljl')
          stateguess.ljl = 0;
        end
        z0(obj.ljl_inds) = stateguess.ljl;
      end
            
      if obj.nC > 0
        %gamma_inds = obj.l_inds(obj.nD+2:obj.nD+2:end, 1);
        lambda_inds = obj.l_inds;          
        if ~isempty(obj.nonlincompl_slack_inds)
          z0(obj.nonlincompl_slack_inds) = obj.nonlincompl_constraint.slack_fun(z0([obj.q_inds;lambda_inds]));
        end
      end
      if ~isempty(obj.jlcompl_slack_inds)
        z0(obj.jlcompl_slack_inds) = obj.jlcompl_constraint.slack_fun(z0([obj.q_inds;obj.ljl_inds]));
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

    function obj = addInputCost(obj,input_cost_function)
      %state_cost.grad_method = 'numerical';
      obj = obj.addCost(input_cost_function,{obj.u_inds});
    end
    
  end

end