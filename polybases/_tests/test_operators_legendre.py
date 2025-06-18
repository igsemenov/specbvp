# -*- coding: utf-8 -*-
"""Test operators on the Legendre polynomials.
"""
import unittest
from specbvp import polybases

LEGENDRE = polybases.Legendre()


class Suite(unittest.TestCase):

    CASES = {}
    OPERATOR = None

    TOL = 1e-11

    def test_suite(self):

        for key in self.CASES:
            self._validate_operator(key)

    def _validate_operator(self, key):

        opr = self.OPERATOR

        point, index, value = self.CASES.get(key)

        opr.setnodes(point)
        opr.setpolys(index)

        self.assertAlmostEqual(
            opr.asdict().get(index), value, msg=key, delta=self.TOL
        )


class TestPolys(Suite):

    OPERATOR = LEGENDRE.polys()

    CASES = {
        'A': (0., 0, 1.),
        'B': (1., 1, 1.),
        'C': (0., 2, -0.5),
        'D': (0., 6, -0.3125)
    }


class TestDerivsOne(Suite):

    OPERATOR = LEGENDRE.derivs(order=1)

    CASES = {
        'A': (1., 0, 0.),
        'B': (0., 1, 1.),
        'C': (0., 3, -0.5*3),
        'D': (0., 7, -0.3125*7)
    }

    PROPS = [
        'DERIV[Pm, 1](0) = m*P_{m-1}(0)'
    ]


class TestDerivsTwo(Suite):

    OPERATOR = LEGENDRE.derivs(order=2)

    CASES = {
        'A': (1., 0, 0.),  # const = 0
        'B': (1., 1, 0.),  # const = 0
        'C': (1., 2, 3.),  # const = 3
        'D': (1., 3, 15),  # 15*x,
        'E': (0., 6, 6*7*0.3125)
    }

    PROPS = [
        'DERIV[Pm, 2](0) = - m*(m+1) * Pm(0)'
    ]


class TestIntegP0Pm(Suite):

    OPERATOR = polybases.legendre.IntegP0Pm()

    CASES = {
        'A': (-1, 0, -2.),
        'B': (0, 1, -0.5),
        'C': (-1, 2, 0),
        'D': (0, 3, 0.125),
        'E': (1., 3, 0.),  # Primint is 0 at x=1.
        'F': (1., 0, 0.),  # Primint is 0 at x=1.
    }


class TestIntegP1Pm(Suite):

    OPERATOR = polybases.legendre.IntegP1Pm()

    CASES = {
        'A': (0., 0, -1./2.),
        'B': (0., 1, -1./3.),
        'C': (0., 2, -1./8.),
        'D': (1., 3, 0.),  # Primint is 0 at x=1.
        'E': (1., 0, 0.),  # Primint is 0 at x=1.
    }


class TestIntegP0PmXB(Suite):

    OPERATOR = LEGENDRE.integxb(weighted=False)

    CASES = {
        'A': (-1, 0, 2.),
        'B': (0, 1, 0.5),
        'C': (0.4, 6, 0.0219408)
    }


class TestIntegP0PmAX(Suite):

    OPERATOR = LEGENDRE.integax(weighted=False)

    CASES = {
        'A': (1., 0, 2.),
        'B': (0., 1, -0.5),
        'C': (0., 5, -1./16.),
        'D': (0., 6, 0)
    }


class TestIntegP1PMXB(Suite):

    OPERATOR = LEGENDRE.integxb(weighted=True)

    CASES = {
        'A': (0., 0, 1./2.),
        'B': (0., 1, 1./3.),
        'C': (0., 2, 1./8.),
        'D': (-0.5, 3, -0.046875),
        'E': (0., 6, 1./128.)
    }


class TestIntegP1PMAX(Suite):

    OPERATOR = LEGENDRE.integax(weighted=True)

    CASES = {
        'A': (0., 0, -1./2.),
        'B': (0., 1, 1./3.),
        'C': (0., 5, 0.),
        'D': (0., 4, 1./48.),
        'E': (0., 6, -1./128.)
    }


if __name__ == '__main__':
    unittest.main()
