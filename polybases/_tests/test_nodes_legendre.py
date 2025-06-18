# -*- coding: utf-8 -*-
"""Test the quadrature rules.
"""
import math
import unittest
import numpy as np
from specbvp import polybases


class GaussRef:
    """Base for the reference Gauss quadratures.
    """

    SOURCE = 'mathworld.wolfram.com/Legendre-GaussGaussrature.html'

    NODES = []
    WEIGHTS = []

    def nodes(self):
        return np.array(self.NODES)

    def weights(self):
        return np.array(self.WEIGHTS)


class GaussTwo(GaussRef):

    NODE = math.sqrt(3.0)/3.0

    NODES = [
        -NODE, NODE
    ]

    WEIGHTS = [
        1., 1.
    ]


class GaussThree(GaussRef):

    NODE = math.sqrt(15.)/5.

    NODES = [
        -NODE, 0., NODE
    ]

    WEIGHTS = [
        5./9., 8./9., 5./9.
    ]


class GaussFour(GaussRef):

    ALPHA = 70.*math.sqrt(30.)

    NODEA = math.sqrt(525. - ALPHA)/35.
    NODEB = math.sqrt(525. + ALPHA)/35.

    WEIGHTA = (18. + math.sqrt(30.))/36.
    WEIGHTB = (18. - math.sqrt(30.))/36.

    NODES = [
        -NODEB, -NODEA, NODEA, NODEB
    ]

    WEIGHTS = [
        WEIGHTB, WEIGHTA, WEIGHTA, WEIGHTB
    ]


class TestGauss(unittest.TestCase):

    TOL = 1e-12

    REFQUADS = {
        2: GaussTwo(),
        3: GaussThree(),
        4: GaussFour()
    }

    QUAD = polybases.Legendre().nodes().get('gauss')

    def assert_equal(self, first, second):

        res = first - second
        err = np.amax(np.fabs(res))

        assert err < self.TOL

    def _validate_quad(self, number):

        quad = self.QUAD.setnum(number)
        quadref = self.REFQUADS.get(number)

        self.assert_equal(
            quad.nodes, quadref.nodes()
        )

        self.assert_equal(
            quad.weights, quadref.weights()
        )

    def test_quad_two(self):
        self._validate_quad(number=2)

    def test_quad_three(self):
        self._validate_quad(number=3)

    def test_quad_four(self):
        self._validate_quad(number=4)


if __name__ == '__main__':
    unittest.main()
