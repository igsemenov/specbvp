# -*- coding: utf-8 -*-
"""Base classes for Newton solvers of non-linear problems.
"""

from abc import ABC, abstractmethod


class IterNewton(ABC):
    """Base class for Newton iterators.
    """

    MAXITER = 20

    def __init__(self):

        self._sol = None
        self._tol = None
        self._count = None
        self._guess = None
        self._maxiter = None
        self._converge = None

    def get_sol(self):
        return self._sol

    def set_config(self, guess, tol, maxiter):
        self.set_tol(tol)
        self.set_guess(guess)
        self.set_maxiter(maxiter)

    def set_tol(self, val):
        self._tol = val

    def set_guess(self, val):
        self._guess = val

    def set_maxiter(self, val):
        self._maxiter = val or self.MAXITER

    def is_converged(self):
        return self._converge

    def __iter__(self):
        self.set_init_state()
        return self

    def set_init_state(self):
        self.set_sol_to_guess()
        self.set_count_to_null()
        self.reset_converge_flag()

    def set_sol_to_guess(self):
        self._sol = self._guess

    def set_count_to_null(self):
        self._count = 0

    def reset_converge_flag(self):
        self._converge = False

    def __next__(self):

        if self._converge is True:
            raise StopIteration
        if self._count == self._maxiter:
            raise StopIteration

        arg = self._sol
        new = self.computenext(arg)

        self._sol = new
        self._count += 1

        res = self.get_residual(new, arg)

        if res < self._tol:
            self._converge = True

        return res

    def computenext(self, arg):
        func = self.get_func(arg)
        deriv = self.get_deriv(arg)
        return self.updatesol(arg, func, deriv)

    def updatesol(self, arg, func, deriv):
        """Updates the argument and returns the result.
        """
        return arg - (func/deriv)

    def get_residual(self, first, second):
        return self.get_norm(first-second)

    @abstractmethod
    def get_norm(self, _):
        """Computes the norm of the argument.
        """

    @abstractmethod
    def get_func(self, _):
        """Computes the function at the argument.
        """

    @abstractmethod
    def get_deriv(self, _):
        """Computes the jacobian at the argument.
        """


class SolverNewton(IterNewton):
    """Base class for Newton solvers.
    """

    def __init__(self):
        super().__init__()
        self.result = None
        self.history = None
        self.converge = None

    def compute(self, guess, tol, maxiter=None):

        self.set_config(guess, tol, maxiter)

        self.run_iterator()
        self.set_converge()
        self.get_results()

        return self

    def run_iterator(self):
        self.history = list(iter(self))

    def get_results(self):
        self.result = self.get_sol()

    def set_converge(self):
        self.converge = self.is_converged()
