# -*- coding: utf-8 -*-
"""Base class for triplet recurrence relations.
"""

from abc import ABC, abstractmethod


class RecurrTriplet(ABC):
    """Base class for triplet recurrence relations.
    """

    def getsequence(self, maxindex) -> list:
        """Computes the recurrence members from 0 to maxindex (>=0).
        """

        startseq = self.genstartseq()

        return list(
            self.runrecurr(startseq, maxindex)
        )

    def runrecurr(self, startseq, maxindex):
        """Generates the recurrence members from 0 to maxindex (>=0).
        """

        if maxindex < len(startseq):
            for item in startseq[0:maxindex+1]:
                yield item
            return

        prev = startseq[-2]
        curr = startseq[-1]

        startindex = len(startseq)-1

        for item in startseq:
            yield item

        for index in range(startindex, maxindex):

            nexter = self.computenext(prev, curr, index)
            yield nexter

            prev = curr
            curr = nexter

    @abstractmethod
    def computenext(self, prev, curr, index):
        """Computes the next member of the recurrence. 
        """

    @abstractmethod
    def genstartseq(self) -> list:
        """Generates the leading recurrence members.
        """
