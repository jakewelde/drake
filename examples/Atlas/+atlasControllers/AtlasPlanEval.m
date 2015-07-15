classdef AtlasPlanEval < PlanEval
% A PlanEval which includes some Atlas-specific tweaks. Specifically, it
% calls its plans' getQPControllerInput method with an additional argument
% (contact_force_detected)

  properties (Access = protected)
    robot
    lcmgl;
  end

  methods
    function obj = AtlasPlanEval(r, varargin)
      obj = obj@PlanEval(varargin{:});
      obj.robot = r;
      obj.lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'atlas');
    end

    function qp_input = getQPControllerInput(obj, t, x, contact_force_detected)
      if nargin < 4
        contact_force_detected = zeros(obj.robot.getNumBodies(), 1);
        contact_force_detected([obj.robot.foot_body_id.right, obj.robot.foot_body_id.left]) = getFootContacts(obj.robot, x(1:obj.robot.getNumPositions));
      end
      plan = obj.getCurrentPlan(t, x);
      qp_input = plan.getQPControllerInput(t, x, contact_force_detected);
    end
  end
end
