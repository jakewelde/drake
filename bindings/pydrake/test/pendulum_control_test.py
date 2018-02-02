from __future__ import print_function
import math
import numpy as np
import unittest


from pydrake.systems.framework import VectorSystem
from pydrake.systems.analysis import Simulator
from pydrake.systems.framework import DiagramBuilder
from pydrake.systems.primitives import SignalLogger

# Define a system to calculate the continuous dynamics
# of a damped pendulum, with mass m, length l,
# gravity g, and damping b.
class DampedPendulumPlant(VectorSystem):
    def __init__(self, m = 3, l = 1, g = 10, b = 2):
        VectorSystem.__init__(self,
            1,                           # One input (torque at shoulder).
            2)                           # One output (theta, dtheta)
        self._DeclareContinuousState(2)  # Two state variables (theta, dtheta).
        self.m = m
        self.l = l
        self.g = g
        self.b = b

    # This method calculates the time derivative of the state,
    # which allows the system to be simulated forward in time.
    # In this case, it implements the continuous time dynamics
    # of a damped pendulum.
    def _DoCalcVectorTimeDerivatives(self, context, u, x, xdot):
        theta = x[0]
        theta_dot = x[1]

        theta_ddot = (  self.m * self.g * self.l * math.sin(theta) \
                      - self.b*theta_dot \
                      + u ) / (self.m * self.l**2)
        xdot[0] = theta_dot
        xdot[1] = theta_ddot

    # This method calculates the output of the system
    # (i.e. those things that are visible downstream of
    # this system) from the state. In this case, it
    # copies out the full state.
    def _DoCalcVectorOutput(self, context, u, x, y):
        y[:] = x

    def _DoHasDirectFeedthrough(self, input_port, output_port):
        return False

class PendulumController(VectorSystem):
    ''' System to control the pendulum. Must be handed
    a function with signature:
        u = f(theta, theta_dot)
    that computes control inputs for the pendulum. '''

    def __init__(self, feedback_rule):
        VectorSystem.__init__(self,
            2,                           # Two inputs (theta, dtheta).
            1)                           # One output (torque at shoulder).
        self.feedback_rule = feedback_rule

    # This method calculates the output of the system from the
    # input by applying the supplied feedback rule.
    def _DoCalcVectorOutput(self, context, u, x, y):
        y[:] = self.feedback_rule(u[0], u[1])



class TestCustom(unittest.TestCase):

    def test_pendulum_control_simulation(x0 = [0.9, 0.0], duration = 10):
        # Create a simple block diagram containing the plant in feedback
        # with the controller.
        pendulum_plant = DampedPendulumPlant()

        # And a controller that returns zero.
        def feedback_rule(theta, dtheta):
            return 0.
        pendulum_controller = PendulumController(feedback_rule)

        builder = DiagramBuilder()
        plant = builder.AddSystem(pendulum_plant)
        controller = builder.AddSystem(pendulum_controller)
        builder.Connect(plant.get_output_port(0), controller.get_input_port(0))
        builder.Connect(controller.get_output_port(0), plant.get_input_port(0))

        # Create a logger to capture the simulation of our plant
        # (We tell the logger to expect a 2-variable input,
        # and hook it up to the pendulum plant's 2-variable output.)
        logger = builder.AddSystem(SignalLogger(2))
        builder.Connect(plant.get_output_port(0), logger.get_input_port(0))

        diagram = builder.Build()

        # Create the simulator.
        simulator = Simulator(diagram)

        # Set the initial conditions for the simulation.
        state = simulator.get_mutable_context().get_mutable_state()\
                         .get_mutable_continuous_state().get_mutable_vector()
        state.SetFromVector(x0)

        # Simulate for the requested duration.
        simulator.StepTo(duration)

if __name__ == '__main__':
    unittest.main()
