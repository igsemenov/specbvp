# -*- coding: utf-8 -*-
"""Test operators on the Chebyschev polynomials.
"""
import math
import unittest
from specbvp import polybases


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

    OPERATOR = polybases.Chebyshev().polys()

    CASES = {
        'A': (0., 0, 1.),
        'B': (1., 1, 1.),
        'C': (
            math.cos(0.5*math.pi/3), 3, 0.
        ),
        'D': (
            math.cos(0.5), 4, math.cos(4*0.5)
        ),
        'E': (0.0, 6, -1.),
        'F': (0.5, 6, +1.)
    }

    PROPS = [
        'Tn[cos(pi/2*n)] = 0',
        'Tn[cos(t)] = cos(n*t)'
    ]


class TestDerivs(Suite):

    OPERATOR = polybases.Chebyshev().derivs(order=1)

    CASES = {
        'A': (1., 0, 0.),
        'B': (0., 1, 1.),
        'C': (
            math.cos(0.5), 4, 4.*math.sin(4*0.5)/math.sin(0.5)
        ),
        'D': (0., 6, 0.),
        'E': (1., 6, 36.)
    }

    PROPS = [
        '(d/dx)T_{n} = n*U_{n-1}',
        'U_{n-1}[cos(t)] = sin(n*t)/sin(t)'
    ]


class TestIntegT0Tn(Suite):

    OPERATOR = polybases.chebyshev.IntegT0Tn()

    CASES = {
        'A': (1., 0, 0.),  # Primint is 0 at x = 1.
        'B': (1., 3, 0.),  # Primint is 0 at x = 1.
        'C': (0.5, 0, -1./2.),
        'D': (0.5, 1, -3./8.),
        'E': (0., 5, -1./6.),
        'F': (0., 6, 1./35.),
        'G': (0., 7, 1./6.),
        'H': (0., 8, 1./63.)
    }


class TestIntegT1Tn(Suite):

    OPERATOR = polybases.chebyshev.IntegT1Tn()

    CASES = {
        'A': (1., 0, 0.),  # Primint is 0 at x = 1.
        'B': (1., 3, 0.),  # Primint is 0 at x = 1.
        'C': (0.5, 0, -3./8.),
        'D': (0., 6, 0.),
        'E': (
            0., 7, (1/2)*(1/35)+(1/2)*(1/63)
        )
    }


class TestIntegT0TnAX(Suite):

    OPERATOR = polybases.Chebyshev().integax(weighted=False)

    CASES = {
        'A': (1., 0, 2.),
        'B': (0., 1, -0.5),
        'C': (0.5, 3, 0.1875),
        'D': (0.5, 4, 0.15)
    }


class TestIntegT0TnXB(Suite):

    OPERATOR = polybases.Chebyshev().integxb(weighted=False)

    CASES = {
        'A': (-1., 0, 2.),
        'B': (0., 1, 0.5),
        'C': (-0.5, 3, -0.1875),
        'D': (-0.5, 4, 0.15)
    }


class TestIntegT1TnAX(Suite):

    OPERATOR = polybases.Chebyshev().integax(weighted=True)

    CASES = {
        'A': (0.5, 0, -3./8.),
        'B': (0.5, 1, +3./8.),
        'C': (0.5, 3, -0.3),
        'D': (0.5, 4, 0.1875)
    }


class TestIntegT1TnXB(Suite):

    OPERATOR = polybases.Chebyshev().integxb(weighted=True)

    CASES = {
        'A': (0.5, 0, +3./8.),
        'B': (0.5, 1, +7./24.),
        'C': (0.5, 3, -0.1),
        'D': (0.5, 4, -0.1875)
    }


if __name__ == '__main__':
    unittest.main()
