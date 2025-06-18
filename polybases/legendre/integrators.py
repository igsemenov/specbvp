# -*- coding: utf-8 -*-
"""Integrate Pm(x) and x*Pm(x) over [x, 1] and [-1, x].
"""

from abc import abstractmethod
import itertools as itr
from . import funcsderivs
from ..abcpolys import PolyOpr

__all__ = [
    'IntegP0Pm', 'IntegP1Pm', 'IntegAX', 'IntegXB'
]


class Primint(PolyOpr):
    """Base class for primitive integrals.
    """

    def getoutputs(self, nodes, maxindex) -> list:
        """Primitive integrals from 0 to maxindex (>=0).
        """

        prim_zero = self.prim_index_zero(nodes)
        prim_from_one = self.prim_index_from_one(nodes, maxindex)

        return self.merge_to_maxindex(
            prim_zero, prim_from_one, maxindex
        )

    def merge_to_maxindex(self, atzero, fromone, maxindex):

        if maxindex == 0:
            return [atzero]

        return [
            atzero, *fromone
        ]

    @abstractmethod
    def prim_index_zero(self, nodes):
        """Primitive integral for m=0.
        """

    @abstractmethod
    def prim_index_from_one(self, nodes, maxindex) -> list:
        """Primitive integral for m >= 1.
        """

    def set_triplets(self, funcs) -> zip:
        """Returns triplets (F_{m-1}, F_{m+1}, m) for m > 1.
        """

        size = len(funcs)

        prev = itr.islice(funcs, 0, size-2)
        coming = itr.islice(funcs, 2, size-0)

        return zip(
            prev, coming, range(1, size-1)
        )


class IntegP0Pm(Primint):
    """Primitive integral of Pm(x) normalized to be 0 at x=1.

    For m = 0:

        INTEGRAL[Pm, x] = x - 1

    For m > 0:

        INTEGRAL[Pm, x] = [P_{m+1} - P_{m-1}]/(2*m+1)

    SOURCE: en.wikipedia.org/wiki/Legendre_polynomials
    """

    POLYS = funcsderivs.Polys()

    def prim_index_zero(self, nodes):
        return nodes - 1.

    def prim_index_from_one(self, nodes, maxindex) -> list:

        if maxindex < 1:
            return []

        polys = self.getpolys(nodes, maxindex+1)
        return self.getprims(polys)

    def getprims(self, polys):

        args = self.set_triplets(polys)

        return list(
            itr.starmap(self.get_integral, args)
        )

    def get_integral(self, prev, coming, count):
        return (coming - prev)/(2*count+1)

    def getpolys(self, nodes, maxindex) -> list:
        return self.POLYS.getoutputs(nodes, maxindex)


class IntegP1Pm(Primint):
    """Primitive integral of x*Pm(x) normalized to be 0 at x=1.

    For m = 0:

        INTEGRAL[x*Pm, x] = (x*x-1)/2

    For m > 0:

        INTEGRAL[x*Pm, x] = (1-ALFA)*BASE_{m+1} + ALFA*BASE_{m-1}

    where:

        BASE_{m} = INTEGRAL[Pm, x]

    with:

        ALFA = m/(2*m+1)

    SOURCE: en.wikipedia.org/wiki/Legendre_polynomials
    REMARK: Derived from the recurrence relation for polynomials.
    """

    def prim_index_zero(self, nodes):
        return 0.5*(nodes*nodes - 1.)

    def prim_index_from_one(self, nodes, maxindex) -> list:

        if maxindex < 1:
            return []

        bases = self.getbases(nodes, maxindex+1)
        return self.getprims(bases)

    def getprims(self, bases):

        args = self.set_triplets(bases)

        return list(
            itr.starmap(self.get_integral, args)
        )

    def get_integral(self, prev, coming, count):
        alfa = self.get_alfa(count)
        return (1.-alfa)*coming + alfa*prev

    def get_alfa(self, count):
        return count/(2*count+1)

    def getbases(self, nodes, maxindex):
        return IntegP0Pm().getoutputs(nodes, maxindex)


class IntegXB(PolyOpr):
    """Base class for integrators over [x, 1].

        INTEGRAL[*, x, 1] = - PRIMINTEG[*](x)

    where PRIMINTEG is anchored as

        PRIMINTEG[*](1) = 0

    """

    PRIMINTEG = None  # Primitive integral that is 0 at x=1.

    def getoutputs(self, nodes, maxindex) -> list:
        outs = self.get_priminteg(nodes, maxindex)
        outs = self.make_negative(outs)
        return outs

    def get_priminteg(self, nodes, maxindex):
        return self.PRIMINTEG.getoutputs(nodes, maxindex)

    def make_negative(self, outs):
        return [
            -val for val in outs
        ]


class IntegAX(PolyOpr):
    """Base class for integrators over [-1, x].

        INTEGRAL[*, -1, x] = TOTAL[*] - INTEGRAL[*, x, 1]

    where

        TOTAL = INTEGRAL[*, -1, 1]

    """

    INTEG_FROM_X = None

    def getoutputs(self, nodes, maxindex) -> list:

        integfromx = self.get_integ_fromx(nodes, maxindex)
        totalinteg = self.get_total_integ(maxindex)

        return self.subtract(
            totalinteg, integfromx
        )

    def get_integ_fromx(self, nodes, maxindex):
        return self.INTEG_FROM_X.getoutputs(nodes, maxindex)

    @abstractmethod
    def get_total_integ(self, maxindex):
        """Integral over [-1, 1] from 0 to maxindex (>=0):
        """

    def subtract(self, itera, iterb):
        return [
            a-b for a, b in zip(itera, iterb)
        ]


class IntegP0PmXB(IntegXB):
    """Integral of Pm(x) over [x, 1].
    """

    PRIMINTEG = IntegP0Pm()


class IntegP0PmAX(IntegAX):
    """Integral of Pm(x) over [-1, x].
    """

    INTEG_FROM_X = IntegP0PmXB()

    def get_total_integ(self, maxindex):
        return [2.] + [0.]*maxindex


class IntegP1PmXB(IntegXB):
    """Integral of x*Pm(x) over [x, 1].
    """

    PRIMINTEG = IntegP1Pm()


class IntegP1PmAX(IntegAX):
    """Integral of x*Pm(x) over [-1, x].
    """

    INTEG_FROM_X = IntegP1PmXB()

    def get_total_integ(self, maxindex):

        if maxindex == 0:
            return [0.]

        if maxindex == 1:
            return [0., 2./3.]

        return [0., 2./3.] + [0.]*(maxindex-1)
