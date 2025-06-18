# -*- coding: utf-8 -*-
"""Integrate Tn(x) and x*Tn(x) over [x, 1] and [-1, x].
"""

from .. import legendre
from . import funcsderivs

__all__ = [
    'IntegT0Tn', 'IntegT1Tn'
]


class IntegT0Tn(legendre.IntegP0Pm):
    """Indefinite integral of Tn(x) normalized to be 0 at x=1.

    For n = 0:

        INTEGRAL[T0, x] = x - 1

    For n = 1:

        INTEGRAL[T1, x] = (x*x - 1)/2

    For n > 1:

        INTEGRAL[Tn, x] = ALFA*T_{n+1} - BETA*T_{n-1} + BIAS

    with:

        ALFA = 1/[2*(n+1)]
        BETA = 1/[2*(n-1)]

        BIAS = 1/[(n+1)*(n-1)]

    SOURCE: en.wikipedia.org/wiki/Chebyshev_polynomials
    """

    POLYS = funcsderivs.Polys()

    def prim_index_zero(self, nodes):
        return nodes - 1.

    def get_integral(self, prev, coming, count):

        if count == 1:
            return 0.25*(coming-prev)

        alfa = self.get_alfa(count)
        beta = self.get_beta(count)
        bias = self.get_bias(count)

        return alfa*coming - beta*prev + bias

    def get_alfa(self, count):
        return 0.5/(count+1.)

    def get_beta(self, count):
        return 0.5/(count-1.)

    def get_bias(self, count):
        return 1.0/(count*count-1.)


class IntegT1Tn(legendre.IntegP1Pm):
    """Indefinite integral of x*Tn(x) normalized to be 0 at x=1.

    For n = 0:

        INTEGRAL[x*Tn, x] = (x*x-1)/2

    For n > 0:

        INTEGRAL[x*Tn, x] = (1/2)*BASE_{n+1} + (1/2)*BASE_{n-1}

    where:

        BASE_{n} = INTEGRAL[Tn, x]

    SOURCE: en.wikipedia.org/wiki/Legendre_polynomials
    REMARK: Derived from the recurrence relation for polynomials.
    """

    def get_integral(self, prev, coming, _):
        return 0.5*prev + 0.5*coming

    def getbases(self, nodes, maxindex):
        return IntegT0Tn().getoutputs(nodes, maxindex)


class IntegT0TnXB(legendre.IntegXB):
    """Integral of Tn(x) over [x, 1].
    """

    PRIMINTEG = IntegT0Tn()


class IntegT0TnAX(legendre.IntegAX):
    """Integral of Tn(x) over [-1, x].
    """

    INTEG_FROM_X = IntegT0TnXB()

    def get_total_integ(self, maxindex):

        if maxindex == 0:
            return [2.]

        if maxindex == 1:
            return [2., 0.]

        def integ(count):
            return (1+(-1)**count)/(1.-count*count)

        return [
            2., 0., *map(integ, range(2, maxindex+1))
        ]


class IntegT1TnXB(legendre.IntegXB):
    """Integral of x*Tn(x) over [x, 1].
    """

    PRIMINTEG = IntegT1Tn()


class IntegT1TnAX(legendre.IntegAX):
    """Integral of x*Tn(x) over [-1, x].
    """

    INTEG_FROM_X = IntegT1TnXB()

    def get_total_integ(self, maxindex):

        if maxindex == 0:
            return [0.]

        if maxindex == 1:
            return [0., 2./3.]

        if maxindex == 2:
            return [0., 2./3., 0.]

        def base(count):
            return (1+(-1)**count)/(1.-count*count)

        def integ(count):
            return (1/2)*base(count-1) + (1/2)*base(count+1)

        return [
            0., 2./3., 0., *map(integ, range(3, maxindex+1))
        ]
