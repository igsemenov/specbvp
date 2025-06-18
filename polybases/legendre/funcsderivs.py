# -*- coding: utf-8 -*-
"""Computes polynomials and their derivatives.
"""

import math
from ..abcpolys import PolyOpr
from ..utils import RecurrTriplet

__all__ = [
    'Polys', 'Derivs'
]


class Recurr(RecurrTriplet):
    """Recurrence for Legendre polynomials and their derivatives.

        NEXT = alfa * x * CURRENT - beta * PREVIOUS

    where:

        alfa = (2*n+1)/(n-m+1)
        beta = (n+m)/(n-m+1)

    with:

        n — Index of CURRENT.
        m — Order of a derivative (from 0).

    SOURCE: de.wikipedia.org/wiki/Legendre-Polynom
    """

    def __init__(self):
        super().__init__()
        self.nodes = None
        self.order = None

    def getoutputs(self, nodes, maxindex):

        _ = self.setnodes(nodes)
        _ = self.getsequence(maxindex)

        return _

    def setnodes(self, nodes):
        self.nodes = nodes
        return self

    def computenext(self, prev, curr, index):

        nodes = self.nodes
        order = self.order

        alfa = self.get_alfa(index, order)
        beta = self.get_beta(index, order)

        return alfa*(nodes*curr) - beta*prev

    def get_alfa(self, index, order):
        return (2.*index+1.)/(index-order+1.)

    def get_beta(self, index, order):
        return (index+order)/(index-order+1.)

    def genstartseq(self):
        pass


class Polys(Recurr, PolyOpr):
    """Computes the Legendre polynomials.
    """

    def __init__(self):
        super().__init__()
        self.order = 0

    def genstartseq(self):

        nodes = self.nodes

        return [
            self.legendre_zero(nodes),
            self.legendre_one(nodes)
        ]

    def legendre_zero(self, nodes):
        return nodes*0. + 1.

    def legendre_one(self, nodes):
        return nodes


class Derivs(Recurr, PolyOpr):
    """Computes the derivatives of Legendre polynomials.
    """

    def __init__(self, order):
        super().__init__()
        self.order = order

    def genstartseq(self):

        nodes = self.nodes

        derivs_below_order = [
            self.derivs_below_order(nodes) for _ in range(self.order)
        ]

        return [
            *derivs_below_order, self.deriv_at_order(nodes)
        ]

    def derivs_below_order(self, nodes):
        return nodes*0.

    def deriv_at_order(self, nodes):

        order = self.order

        _ = math.factorial(order)*math.pow(2., order)
        deriv = math.factorial(2*order)/_

        return nodes*0. + deriv
