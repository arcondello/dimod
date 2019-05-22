# Copyright 2019 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
# =============================================================================
from __future__ import division

try:
    import collections.abc as abc
except ImportError:
    import collections as abc

from copy import deepcopy
from functools import wraps
from itertools import chain

import numpy as np

from dimod.bqm.cyutils import coo_sort
from dimod.decorators import vartype_argument
from dimod.exceptions import WriteableError
from dimod.variables import Variables
from dimod.vartypes import BINARY, SPIN

__all__ = ['CooBinaryQuadraticModel',
           'CooBQM',
           ]


class DenseLinear(abc.Mapping):
    def __init__(self, bqm):
        self.bqm = bqm

    def __getitem__(self, v):
        return self.bqm.ldata[self.bqm.variables.index[v]]

    def __setitem__(self, v, bias):

        # developer note: The user can manually set bqm.ldata.flags.writeable
        # and that would make ths function work even if bqm.is_writeable is
        # False. However, I think it makes sense for .linear to be consistent
        # with the rest of the bqm rather than the underlying data array only.
        # See also CooQuadratic.__setitem__
        if not self.bqm.is_writeable:
            msg = 'Cannot set linear bias when {}.is_writeable is False'
            raise WriteableError(msg.format(type(self.bqm).__name__))

        self.bqm.ldata[self.bqm.variables.index[v]] = bias

    def __iter__(self):
        return iter(self.bqm.variables)

    def __len__(self):
        return len(self.bqm.variables)

    def __repr__(self):
        # todo: not cast to dict first
        return str(dict(self))


class CooQuadratic(abc.Mapping):
    def __init__(self, bqm):
        self.bqm = bqm

    def _get_index(self, interaction):
        # get either a boolean mask or a slice for a particular interaction
        bqm = self.bqm
        irow = bqm.irow
        icol = bqm.icol
        variables = bqm.variables

        u, v = interaction

        iu = variables.index[u]
        iv = variables.index[v]

        if bqm.is_sorted:
            # O(log(|E|))
            if iu > iv:
                iu, iv = iv, iu

            # we want the window in irow. There are more clever ways to get
            # the right side (going back up the binary search tree) but this
            # is good enough for now
            row_left = np.searchsorted(irow, iu, side='left')
            row_right = row_left + np.searchsorted(irow[row_left:], iu,
                                                   side='right')

            # there may be duplicates so again we need to check for a window
            col_left = row_left + np.searchsorted(icol[row_left:row_right], iv,
                                                  side='left')

            # if we know it is aggregated we could skip this step
            col_right = col_left + np.searchsorted(irow[col_left:row_right], iv,
                                                   side='right')

            if irow[col_left] == iu and icol[col_left] == iv:
                return slice(col_left, col_right)
        else:
            # O(|E|)
            # if we cythonized it we could avoid making the mask which would
            # save on memory. Creating the mask is faster than iterating in
            # python with a list
            mask = ((irow == iu) & (icol == iv)) ^ ((irow == iv) & (icol == iu))
            if np.any(mask):
                return mask

        raise KeyError

    def __getitem__(self, interaction):
        return self.bqm.qdata[self._get_index(interaction)].sum()

    def __setitem__(self, interaction, bias):
        # developer note: See note in DenseLinear.__setitem__
        if not self.bqm.is_writeable:
            msg = 'Cannot set linear bias when {}.is_writeable is False'

        index = self._get_index(interaction)

        qdata = self.bqm.qdata

        # because there might be duplicates, we set them all to zero before
        # setting the first equal to bias
        qdata[index] = 0
        qdata[index][0] = bias  # _get_index guarantees this exists

    def __iter__(self):
        bqm = self.bqm
        variables = bqm.variables
        for iu, iv in zip(bqm.irow, bqm.icol):
            yield variables[iu], variables[iv]

    def __len__(self):
        return len(self.bqm.qdata)

    def __repr__(self):
        # todo: not cast to dict first
        return str(dict(self.items()))

    def items(self):
        return CooQuadraticItemsView(self)  # much faster to iterate directly


class CooQuadraticItemsView(abc.ItemsView):
    def __iter__(self):
        bqm = self._mapping.bqm
        variables = bqm.variables
        for iu, iv, bias in zip(bqm.irow, bqm.icol, bqm.qdata):
            yield (variables[iu], variables[iv]), bias


class CooBinaryQuadraticModel(object):
    """

    Notes:
        The linear biases are stored in "dense" form.


    """
    @vartype_argument('vartype')
    def __init__(self, linear, quadratic, offset, vartype,
                 variables=None,
                 dtype=np.float, index_dtype=np.int,
                 copy=True):

        self.dtype = dtype = np.dtype(dtype)
        self.index_dtype = index_dtype = np.dtype(index_dtype)

        # linear
        self.ldata = ldata = np.array(linear, dtype=dtype, copy=copy)

        # quadratic
        irow, icol, qdata = quadratic
        self.irow = irow = np.array(irow, dtype=index_dtype, copy=copy)
        self.icol = icol = np.array(icol, dtype=index_dtype, copy=copy)
        self.qdata = qdata = np.array(qdata, dtype=dtype, copy=copy)

        # CooBQMs have fixed shape
        irow.flags.writeable = False
        icol.flags.writeable = False

        # todo: check that they are all 1dim, same length etc

        # offset
        self.offset = dtype.type(offset)

        # vartype
        self.vartype = vartype

        # variables
        if not copy and isinstance(variables, Variables):
            self.variables = variables
        elif variables is None:
            self.variables = Variables(range(len(ldata)))
        else:
            self.variables = Variables(variables)

        # views
        self.linear = DenseLinear(self)
        self.quadratic = CooQuadratic(self)

    def __eq__(self, other):
        # todo: performance here can be much improved
        return (self.linear == other.linear and
                self.quadratic == other.quadratic and
                self.offset == other.offset and
                self.vartype == other.vartype)

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return '{}.from_dicts({}, {}, {}, {!r})'.format(type(self).__name__,
                                                        self.linear,
                                                        self.quadratic,
                                                        self.offset,
                                                        self.vartype.name)

    #
    # Flags
    #

    @property
    def is_sorted(self):
        """True if irow, icol are sorted lexicographically"""

        # because the irow/icol are immutable, we can keep the results cached
        if hasattr(self, '_is_sorted'):
            return self._is_sorted

        # ideally we could speed this up with cython, but right now
        # memoryviews don't handle read-only arrays very well (see
        # https://github.com/dask/distributed/issues/1978) so we just do this
        # in numpy
        irow = self.irow
        icol = self.icol

        is_sorted = ((irow <= icol).all() and
                     (irow[:-1] <= irow[1:]).all() and
                     (icol[:-1] <= icol[1:])[irow[:-1] == irow[1:]].all())

        self._is_sorted = is_sorted
        return is_sorted

    @property
    def is_writeable(self):
        return (self.variables.is_writeable or
                self.ldata.flags.writeable or
                self.qdata.flags.writeable)

    def setflags(self, write=None):
        if write is not None:
            self.ldata.flags.writeable = self.qdata.flags.writeable = write
            self.variables.setflags(write=write)

    #
    # Construction
    #

    @classmethod
    def from_dicts(cls, linear, quadratic, offset, vartype,
                   dtype=np.float, index_dtype=np.int):

        # get all of the labels
        variables = Variables(chain(linear, chain(*quadratic)))

        # get the quadratic
        num_interactions = len(quadratic)

        irow = np.empty(num_interactions, dtype=index_dtype)
        icol = np.empty(num_interactions, dtype=index_dtype)
        qdata = np.empty(num_interactions, dtype=dtype)

        # we could probably speed this up with cython
        for idx, ((u, v), bias) in enumerate(quadratic.items()):
            irow[idx] = variables.index[u]
            icol[idx] = variables.index[v]
            qdata[idx] = bias

        # next linear
        ldata = np.fromiter((linear.get(v, 0) for v in variables),
                            count=len(variables), dtype=dtype)

        bqm = cls(ldata, (irow, icol, qdata), offset, vartype,
                  variables=variables,
                  dtype=dtype, index_dtype=index_dtype,
                  copy=False)  # we just made these objects so no need to copy

        return bqm

    @classmethod
    def from_ising(cls, h, J, offset=0,
                   dtype=np.float, index_dtype=np.int):
        return cls.from_dicts(h, J, offset, SPIN,
                              dtype=dtype, index_dtype=index_dtype)

    @classmethod
    def from_qubo(cls, Q, offset=0, dtype=np.float, index_dtype=np.int):
        return cls.from_dicts({}, Q, offset, BINARY,
                              dtype=dtype, index_dtype=index_dtype)

    #
    # Methods
    #

    # developer note: this method returns a new object in order to keep it
    # consistent with SampleSet.aggregate and because unlike .relabel
    # and .sort it changes the size of the bqm.
    def aggregate(self, copy=True):
        """Return a binary quadratic model with any duplicate interactions
        aggregated.

        Note that this will also sort the interaction vectors.

        Args:
            copy (bool, optional, default=True):
                If True a new binary quadratic model is always returned. If
                False a copy is made only when there are duplicate interactions.

        Returns:
            :obj:`.CooBinaryQuadraticModel`

        """
        if self.is_sorted:
            base = self
        else:
            base = self.copy()
            base.sort()
            copy = False  # we've done a copy so no need for another

        irow = base.irow
        icol = base.icol

        # this is a case that can be sped-up with cython to avoid the
        # intermediate objects. Though it is still faster than doing python
        # loops
        unique = np.empty(len(irow), dtype=bool)
        if len(irow):
            unique[0] = True
        unique[1:] = ((irow[1:] != irow[:-1]) | (icol[1:] != icol[:-1]))

        if unique.all():
            qdata = base.qdata
        else:
            irow = irow[unique]
            icol = icol[unique]

            unique_idxs, = np.nonzero(unique)
            qdata = np.add.reduceat(base.qdata, unique_idxs, dtype=base.dtype)

            copy = False  # we've done a copy so no need for another

        return CooBQM(base.ldata,
                      (irow, icol, qdata),
                      base.offset,
                      base.vartype,
                      variables=base.variables,
                      dtype=base.dtype, index_dtype=base.index_dtype,
                      copy=copy)

    @vartype_argument('vartype')
    def change_vartype(self, vartype):
        """Change the binary quadratic model's vartype in-place."""
        if not self.is_writeable:
            msg = "cannot be sorted while {}.is_writeable is set to False"
            raise WriteableError(msg.format(type(self).__name__))

        if self.vartype is BINARY and vartype is SPIN:

            # qx = q(s+1)/2
            self.ldata /= 2
            self.offset += self.ldata.sum()

            # x'Qx = (s'Qs + 1'Qs + s'Q1 + 1'Q1) / 4
            self.qdata /= 4

            np.add.at(self.ldata, self.irow, self.qdata)
            np.add.at(self.ldata, self.icol, self.qdata)

            self.offset += self.qdata.sum()

            self.vartype = SPIN

        elif self.vartype is SPIN and vartype is BINARY:

            # hs = 2hx - h1
            self.offset -= self.ldata.sum()
            self.ldata *= 2

            # s'Js = 4x'Jx - 2*1'Jx - 2*x'J1 + 1'J1
            self.offset += self.qdata.sum()

            # we want qdata *= 4, but to save on making an intermediate object
            # in the next step, we do half of it here
            self.qdata *= 2

            np.subtract.at(self.ldata, self.irow, self.qdata)
            np.subtract.at(self.ldata, self.icol, self.qdata)

            self.qdata *= 2

            self.vartype = BINARY

    def copy(self):
        return deepcopy(self)

    def relabel(self, mapping):
        try:
            self.variables.relabel(mapping)
        except WriteableError as e:
            msg = "cannot be relabelled while {}.is_writeable is set to False"
            e.args = (msg.format(type(self).__name__), e.args[1:])
            raise e

    def sort(self):
        """Sort the interactions in-place.

        Note that sort is not stable.
        """
        if self.is_sorted:
            return  # we are already done

        if not self.is_writeable:
            msg = "cannot be sorted while {}.is_writeable is set to False"
            raise WriteableError(msg.format(type(self).__name__))

        # make the irow/icol writeable
        self.irow.flags.writeable = True
        self.icol.flags.writeable = True

        coo_sort(self.irow, self.icol, self.qdata)  # acts in-place

        # mark them not writeable again
        self.irow.flags.writeable = False
        self.icol.flags.writeable = False

        # we're now sorted
        self._is_sorted = True


CooBQM = CooBinaryQuadraticModel
