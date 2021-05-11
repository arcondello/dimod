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

include "cybqm_template.pyx.pxi"

import numpy as np

BIAS_DTYPE = np.float64
INDEX_DTYPE = np.int64

__all__ = ['cyBQM_float64']

cdef class cyBQM_float64(cyBQM_template):
    pass
