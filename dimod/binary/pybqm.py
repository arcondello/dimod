# Copyright 2021 D-Wave Systems Inc.
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

from __future__ import annotations

import copy

from collections.abc import Collection, Mapping, KeysView, Callable
from functools import reduce
from typing import Any, Dict, Iterator, Optional, Tuple

import numpy as np

try:
    from numpy.typing import ArrayLike, DTypeLike
except ImportError:
    ArrayLike = Any
    DTypeLike = Any

from dimod.sampleset import as_samples
from dimod.typing import Variable, VartypeLike
from dimod.utilities import iter_safe_relabels
from dimod.vartypes import as_vartype, Vartype


class pyBQM:
    """A pure-python BQM implementation for handling arbitrary bias types."""
    def __init__(self, vartype: VartypeLike):
        self._adj: Dict[Variable, Dict[Variable, Any]] = dict()
        self._vartype = as_vartype(vartype)

        self.offset = 0.0

    @property
    def dtype(self) -> np.dtype:
        return np.dtype('O')

    @property
    def variables(self) -> Collection:
        return KeysView(self._adj)

    @property
    def vartype(self):
        return self._vartype

    def add_linear(self, v: Variable, bias: Any):
        self._adj.setdefault(v, dict())
        self._adj[v][v] = self._adj[v].get(v, 0) + bias

    def add_linear_from_array(self, linear: ArrayLike):
        for v, bias in enumerate(np.asarray(linear)):
            self.add_linear(v, bias)

    def add_quadratic(self, u: Variable, v: Variable, bias: Any):
        if u == v:
            raise ValueError(f"{u!r} cannot have an interaction with itself")

        self.add_variable(u)
        self.add_variable(v)

        zero = type(bias)(0)

        self._adj[u][v] = self._adj[v][u] = bias + self._adj[v].get(u, zero)

    def add_quadratic_from_dense(self, quadratic: ArrayLike):
        quadratic = np.asarray(quadratic)
        if quadratic.shape[0] != quadratic.shape[1]:
            raise ValueError("quadratic must be a square matrix")

        num_variables = quadratic.shape[0]

        for u in range(num_variables):
            for v in range(num_variables):
                if u == v:
                    self.add_linear(u, quadratic[u, v])
                elif quadratic[u, v]:
                    self.add_quadratic(u, v, quadratic[u, v])

    def add_variable(self, v: Optional[Variable] = None,
                     bias: Any = 0) -> Variable:
        if v is None:
            # match the behaviour of dimod.Variables
            adj = self._adj

            v = len(adj)

            if v in adj:
                v = 0
                while v in adj:
                    v += 1

        self.add_linear(v, bias)
        return v

    def change_vartype(self, vartype: VartypeLike) -> pyBQM:
        vartype = as_vartype(vartype)

        # in place and we are already correct, so nothing to do
        if self.vartype == vartype:
            return self

        if vartype == Vartype.BINARY:
            lin_mp = 2.
            lin_offset_mp = -1.
            quad_mp = 4.
            lin_quad_mp = -2.
            quad_offset_mp = .5
        elif vartype == Vartype.SPIN:
            lin_mp = .5
            lin_offset_mp = .5
            quad_mp = .25
            lin_quad_mp = .25
            quad_offset_mp = .125
        else:
            raise RuntimeError("unexpected vartype")

        adj = self._adj

        for u, Nu in adj.items():
            lbias = Nu[u]

            self.offset += lin_offset_mp * lbias
            Nu[u] = lin_mp * lbias

            for v, qbias in Nu.items():
                if v == u:
                    continue
                Nu[v] = quad_mp * qbias
                Nu[u] += lin_quad_mp * qbias  # linear
                self.offset += quad_offset_mp * qbias

        self._vartype = vartype

        return self

    def degree(self, v: Variable) -> int:
        try:
            return len(self._adj[v]) - 1
        except KeyError:
            raise ValueError(f"unknown variable {v!r}") from None

    def energies(self, samples_like, dtype: DTypeLike = None):
        samples, labels = as_samples(samples_like)

        ldata, (irow, icol, qdata), offset \
            = self.to_numpy_vectors(variable_order=labels)

        energies = samples.dot(ldata)
        energies += (samples[:, irow]*samples[:, icol]).dot(qdata)
        energies += offset
        return np.asarray(energies, dtype=dtype)

    def get_linear(self, v: Variable) -> Any:
        try:
            return self._adj[v][v]
        except KeyError:
            raise ValueError(f"unknown variable {v!r}") from None

    def get_quadratic(self, u: Variable, v: Variable,
                      default: Optional[Any] = None) -> Any:
        if u == v:
            raise ValueError(f"{u!r} cannot have an interaction with itself")

        try:
            return self._adj[u][v]
        except KeyError:
            if default is None:
                raise ValueError(
                    f"{u!r} and {v!r} have no interaction") from None
            return default

    def iter_neighborhood(self, v: Variable) -> Iterator[Any]:
        try:
            Nv = self._adj[v]
        except KeyError:
            raise ValueError(f"unknown variable {v!r}") from None
        for u, bias in Nv.items():
            if u != v:
                yield u, bias

    def iter_quadratic(self) -> Iterator[Tuple[Variable, Variable, Any]]:
        seen = set()
        for u, Nu in self._adj.items():
            seen.add(u)
            for v, bias in Nu.items():
                if v not in seen:
                    yield u, v, bias

    def num_variables(self) -> int:
        return len(self._adj)

    def num_interactions(self) -> int:
        n = sum(map(len, self._adj.values()), 0)
        n -= self.num_variables()  # subtract the self-loops
        return n // 2

    def reduce_linear(self, function: Callable,
                      initializer: Optional[Any] = None) -> Any:
        gen = (self.get_linear(v) for v in self.variables)
        if initializer is None:
            return reduce(function, gen)
        else:
            return reduce(function, gen, initializer)

    def reduce_neighborhood(self, v: Variable, function: Callable,
                            initializer: Optional[Any] = None) -> Any:
        gen = (b for _, b in self.iter_neighborhood(v))
        if initializer is None:
            return reduce(function, gen)
        else:
            return reduce(function, gen, initializer)

    def reduce_quadratic(self, function: Callable,
                         initializer: Optional[Any] = None) -> Any:
        gen = (b for _, _, b in self.iter_quadratic())
        if initializer is None:
            return reduce(function, gen)
        else:
            return reduce(function, gen, initializer)

    def remove_interaction(self, u: Variable, v: Variable):
        if u == v:
            raise ValueError(f"{u!r} cannot have an interaction with itself")
        try:
            self._adj[u].pop(v)
        except KeyError:
            raise ValueError(f"{u!r} and {v!r} have no interaction") from None

        self._adj[v].pop(u)

    def relabel_variables(self, mapping: Mapping[Variable, Variable]):
        adj = self._adj

        for submap in iter_safe_relabels(mapping, self.variables):
            for old, new in submap.items():
                if old == new:
                    continue

                # replace the linear bias
                adj[new] = {new: adj[old].pop(old)}

                # copy the quadratic biases
                for v in adj[old]:
                    adj[new][v] = adj[v][new] = adj[v].pop(old)

                # remove the old adj for old
                del adj[old]

    def relabel_variables_as_integers(self) -> Mapping[int, Variable]:
        mapping = dict((v, i) for i, v in enumerate(self.variables) if i != v)
        self.relabel_variables(mapping)
        return dict((i, v) for v, i in mapping.items())

    def remove_variable(self, v: Optional[Variable] = None) -> Variable:
        if v is None:
            try:
                v, Nv = self._adj.popitem()
            except KeyError:
                raise ValueError("cannot pop from an empty model") from None
        else:
            try:
                Nv = self._adj.pop(v)
            except KeyError:
                raise ValueError(f"unknown variable {v!r}") from None

        for u in Nv:
            if u != v:
                self._adj[u].pop(v)

        return v

    def resize(self, n: int):
        while n > self.num_variables():
            self.add_variable()
        while n < self.num_variables():
            self.remove_variable()

    def set_linear(self, v: Variable, bias: Any):
        self._adj.setdefault(v, dict())[v] = bias

    def set_quadratic(self, u: Variable, v: Variable, bias: Any):
        if u == v:
            raise ValueError(f"{u!r} cannot have an interaction with itself")
        self.add_variable(u)
        self.add_variable(v)
        self._adj[u][v] = self._adj[v][u] = bias

    def to_numpy_vectors(self, variable_order=None,
                         sort_indices=False, sort_labels=True,
                         return_labels=False):
        num_variables = self.num_variables()
        num_interactions = self.num_interactions()

        if variable_order is None:
            variable_order = list(self._adj)

            if sort_labels:
                try:
                    variable_order.sort()
                except TypeError:
                    # can't sort unlike types in py3
                    pass

        if len(variable_order) != num_variables:
            raise ValueError("variable_order does not match the number of "
                             "variables")

        ldata = np.asarray([self.get_linear(v) for v in variable_order])

        label_to_idx = {v: idx for idx, v in enumerate(variable_order)}
        irow = []
        icol = []
        qdata = []
        for u, v, bias in self.iter_quadratic():
            irow.append(label_to_idx[u])
            icol.append(label_to_idx[v])
            qdata.append(bias)

        irow = np.asarray(irow)
        icol = np.asarray(icol)
        qdata = np.asarray(qdata)

        if sort_indices:
            # row index should be less than col index, this handles
            # upper-triangular vs lower-triangular
            swaps = irow > icol
            if swaps.any():
                # in-place
                irow[swaps], icol[swaps] = icol[swaps], irow[swaps]

            # sort lexigraphically
            order = np.lexsort((irow, icol))
            if not (order == range(len(order))).all():
                # copy
                irow = irow[order]
                icol = icol[order]
                qdata = qdata[order]

        ret = [ldata, (irow, icol, qdata), ldata.dtype.type(self.offset)]

        if return_labels:
            ret.append(variable_order)

        return tuple(ret)

    def update(self, other):
        # in the future when we have vartype views we can improve this
        if other.vartype != self.vartype:
            other = copy.deepcopy(other)
            other.change_vartype(self.vartype)

        for v in other.variables:
            self.add_linear(v, other.get_linear(v))

        for u, v, bias in other.iter_quadratic():
            self.add_quadratic(u, v, bias)

        self.offset += other.offset
