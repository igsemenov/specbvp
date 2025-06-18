# -*- coding: utf-8 -*-
"""Public class for a polynomial basis.
"""

from ..abcpolys import PolyBasis
from . import funcsderivs
from . import integrators
from . import gaussnodes


def apiobj(obj):
    obj.__module__ = 'polybases'
    return obj


def mutedocs(obj):
    for val in vars(obj).values():
        if callable(val):
            val.__doc__ = ''
    return obj


@apiobj
@mutedocs
class Legendre(PolyBasis):
    """Basis formed by the Legendre polynomials.
    """

    def polys(self):
        return funcsderivs.Polys()

    def derivs(self, order):
        return funcsderivs.Derivs(order)

    def integax(self, weighted=False):
        if not weighted:
            return integrators.IntegP0PmAX()
        return integrators.IntegP1PmAX()

    def integxb(self, weighted=False):
        if not weighted:
            return integrators.IntegP0PmXB()
        return integrators.IntegP1PmXB()

    def nodes(self):
        return {
            'gauss': gaussnodes.GaussNodes()
        }
