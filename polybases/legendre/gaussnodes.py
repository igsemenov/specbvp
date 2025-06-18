# -*- coding: utf-8 -*-
"""Class representing the Gauss nodes in `(-1,1)`. 
"""

import numpy as np
from ..abcpolys import NodeSet
from ..utils import findroots
from . import funcsderivs


class GaussNodes(NodeSet):
    """Set of Gauss nodes in [-1, 1].
    """

    TOL = 1e-14
    MAXITER = 10

    def find_nodes(self, number):

        finder = NodesFinder(number)
        guess = finder.get_nodes_guess(number)

        finder.compute(
            guess, tol=self.TOL, maxiter=self.MAXITER
        )

        if finder.converge is True:
            return finder.result

        raise NodesError(
            'Gauss nodes finder failed, no convergence'
        )

    def find_weights(self, nodes):
        return WeightsFinder().compute_weights(nodes)


class NodesFinder(findroots.SolverNewton):
    """Computes the Gauss nodes.

        x = x - POLYN(x)/DERIVN(x)

    where:

        POLYN — Legendre polynomial of the n-th order.
        DERIVN — First derivative of POLYN.

    The solution guess is

        x = COS[pi*(4*k-1)/(4*n+2)]

    with

        k = n, ..., 1 (order matters)

    SOURCE: de.wikipedia.org/wiki/Legendre-Polynom#Nullstellen
    """

    POLYS = funcsderivs.Polys()
    DERIVS = funcsderivs.Derivs(order=1)

    def __init__(self, number):
        super().__init__()
        self.number = number

    def get_norm(self, res):

        maxval = res.max()
        minval = res.min()

        return max(abs(maxval), abs(minval))

    def get_func(self, nodes):

        outs = self.POLYS.getoutputs(
            nodes=nodes, maxindex=self.number
        )

        return outs.pop()

    def get_deriv(self, nodes):

        outs = self.DERIVS.getoutputs(
            nodes=nodes, maxindex=self.number
        )

        return outs.pop()

    def get_nodes_guess(self, number):
        index = np.arange(number, 0, -1)
        return np.cos(np.pi*(4*index-1.)/(4*number+2.))


class WeightsFinder:
    """Computes weights of the Gauss quadrature.

        w = 2/[n*POLYN1(x)*DERIVN(x)]

    where:

        POLYN1 — Legendre polynomial of the order n-1.
        DERIVN — First derivative of the n-th Legendre polynomial.

    with:

        n — Number of the quadrature nodes. 
        x — Nodes of the Gauss-Legendre quadrature.

    SOURCE: mathworld.wolfram.com/Legendre-GaussQuadrature.html
    """

    POLYS = funcsderivs.Polys()
    DERIVS = funcsderivs.Derivs(order=1)

    def compute_weights(self, nodes):

        number = nodes.shape[0]

        polys = self.getpolys(nodes, number-1)
        derivs = self.getderivs(nodes, number-0)

        return 2.0/(number*polys*derivs)

    def getpolys(self, nodes, index):
        outs = self.POLYS.getoutputs(nodes, index)
        return outs.pop()

    def getderivs(self, nodes, index):
        outs = self.DERIVS.getoutputs(nodes, index)
        return outs.pop()


class NodesError(Exception):
    """Raised when the nodes finder does not converge.
    """
