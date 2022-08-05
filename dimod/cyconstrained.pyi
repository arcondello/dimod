# Copyright 2022 D-Wave Systems Inc.
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

import typing

import dimod


class cyConstrainedQuadraticModel:
    def add_constraint_from_model(self,
                                  qm: typing.Union[dimod.binary.cybqm.cyBQM_float32,
                                                   dimod.binary.cybqm.cyBQM_float64,
                                                   dimod.quadratic.cyqm.cyQM_float32,
                                                   dimod.quadratic.cyqm.cyQM_float64],
                                  sense: typing.Union[dimod.sym.Sense, str],
                                  rhs: dimod.typing.Bias = 0,
                                  *,
                                  label: typing.Optional[typing.Hashable] = None,
                                  copy: bool = True) -> typing.Hashable: ...

    def add_variable(self,
                     vartype: dimod.typing.VartypeLike,
                     v: typing.Optional[dimod.typing.Variable],
                     *,
                     lower_bound: typing.Optional[float] = None,
                     upper_bound: typing.Optional[float] = None,
                     ) -> dimod.typing.Variable: ...

    def lower_bound(self, v: dimod.typing.Variable) -> float: ...
    def set_lower_bound(self, v: dimod.typing.Variable, lb: float) -> None: ...

    def set_objective(self,
                      objective: typing.Union[dimod.binary.cybqm.cyBQM_float32,
                                              dimod.binary.cybqm.cyBQM_float64,
                                              dimod.quadratic.cyqm.cyQM_float32,
                                              dimod.quadratic.cyqm.cyQM_float64],
                      ) -> None: ...

    def set_upper_bound(self, v: dimod.typing.Variable, ub: float) -> None: ...
    def upper_bound(self, v: dimod.typing.Variable) -> float: ...
    def vartype(self, v: dimod.typing.Variable) -> dimod.Vartype: ...
