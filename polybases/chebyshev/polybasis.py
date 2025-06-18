# -*- coding: utf-8 -*-
"""Public class for a polynomial basis.
"""

from ..abcpolys import PolyBasis
from . import funcsderivs
from . import integrators


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
class Chebyshev(PolyBasis):
    """Basis formed by the Chebyshev polynomials of the 1st kind.
    """

    def polys(self):
        return funcsderivs.Polys()

    def derivs(self, order):
        if order == 1:
            return funcsderivs.Derivs()
        raise NotImplementedError(
            'currently not implemented for Chebyshev'
        )

    def integax(self, weighted=False):
        if not weighted:
            return integrators.IntegT0TnAX()
        return integrators.IntegT1TnAX()

    def integxb(self, weighted=False):
        if not weighted:
            return integrators.IntegT0TnXB()
        return integrators.IntegT1TnXB()

    def nodes(self):
        raise NotImplementedError(
            'currently not implemented for Chebyshev'
        )
