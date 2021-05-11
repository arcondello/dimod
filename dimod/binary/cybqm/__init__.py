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

import numpy as np

from dimod.binary.cybqm.cybqm_float32 import cyBQM_float32
from dimod.binary.cybqm.cybqm_float64 import cyBQM_float64

__all__ = [
    'cyBQM_float32',
    'cyBQM_float64',
    ]


def cybqm_cls(dtype):
    if np.issubdtype(dtype, np.float32):
        return cyBQM_float32
    elif np.issubdtype(dtype, np.float64):
        return cyBQM_float64
    else:
        raise ValueError("unsupported dtype")
