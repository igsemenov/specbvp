# -*- coding: utf-8 -*-
"""Test Newton solvers for non-linear problems.
"""
import math
import unittest
from specbvp.polybases import utils


class DottieNumber(utils.SolverNewton):
    """Example Newton solver derived from the base one.
    """

    def get_norm(self, arg):
        return abs(arg)

    def get_func(self, arg):
        return math.cos(arg) - arg

    def get_deriv(self, arg):
        return -math.sin(arg) - 1.


class TestSolverNewton(unittest.TestCase):

    DOTTIE = 0.7390851332151607

    def test_compute(self):

        solver = DottieNumber()

        solver.compute(
            guess=0., tol=1e-8
        )

        self.assertAlmostEqual(
            solver.result, self.DOTTIE, delta=1e-8
        )

    def test_converge(self):

        solver = DottieNumber()

        assert solver.compute(0., 1e-6).converge is True
        assert solver.compute(0., 1e-6, maxiter=1).converge is False

    def test_history(self):

        solver = DottieNumber()
        history = solver.compute(guess=0., tol=1e-6, maxiter=3).history

        self.assertAlmostEqual(
            history[0], 1.0, delta=1e-10
        )

        self.assertAlmostEqual(
            history[1], 0.249, delta=1e-3
        )


if __name__ == '__main__':
    unittest.main()
