# -*- coding: utf-8 -*-
"""Computes polynomials and their derivatives.
"""

from abc import abstractmethod
from ..abcpolys import PolyOpr
from ..utils import RecurrTriplet


class Recurr(RecurrTriplet):
    """Recurrence for Chebyschev polynomials of both kinds.

        NEXT = 2 * x * CURRENT - PREVIOUS

    SOURCE: en.wikipedia.org/wiki/Chebyshev_polynomials
    """

    def __init__(self):
        super().__init__()
        self.nodes = None

    def setnodes(self, nodes):
        self.nodes = nodes
        return self

    def computenext(self, prev, curr, _):
        return 2.*(self.nodes*curr) - prev

    def genstartseq(self):

        nodes = self.nodes

        return [
            self.poly_zero(nodes),
            self.poly_one(nodes)
        ]

    @abstractmethod
    def poly_zero(self, nodes):
        pass

    @abstractmethod
    def poly_one(self, nodes):
        pass


class ChebOne(Recurr):
    """Recurrence for polynomials of the 1st kind.

    T0 = 1
    T1 = x
    T2 = ...

    """

    def poly_zero(self, nodes):
        return nodes*0. + 1.

    def poly_one(self, nodes):
        return nodes


class ChebTwo(Recurr):
    """Recurrence for polynomials of the 2nd kind.

    U0 = 1
    U1 = 2*x
    U2 = ...

    """

    def poly_zero(self, nodes):
        return nodes*0. + 1.

    def poly_one(self, nodes):
        return 2.*nodes


class Polys(ChebOne, PolyOpr):
    """Operator for getting basis polynomials — Tn(x).
    """

    def getoutputs(self, nodes, maxindex):

        _ = self.setnodes(nodes)
        _ = self.getsequence(maxindex)

        return _


class Derivs(ChebTwo, PolyOpr):
    """Operator for getting derivatives of the basis polynomials.

        DERIV[Tn, x] = n*U_{n-1}

    SOURCE: en.wikipedia.org/wiki/Chebyshev_polynomials
    """

    def getoutputs(self, nodes, maxindex) -> list:

        if maxindex == 0:
            return [
                nodes*0.
            ]

        return self.getderivs(nodes, maxindex)

    def getupolys(self, nodes, maxindex) -> list:

        _ = self.setnodes(nodes)
        _ = self.getsequence(maxindex)

        return _

    def getderivs(self, nodes, maxindex) -> list:

        upolys = self.getupolys(nodes, maxindex)

        derivs_from_one = [
            index*upolys[index-1] for index in range(1, maxindex+1)
        ]

        return [
           nodes*0., *derivs_from_one
        ]
