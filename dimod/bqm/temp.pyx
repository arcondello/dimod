from dimod.bqm.cppbqm cimport cppAdjMapBQM
from dimod.bqm.adjarraybqm cimport cyAdjArrayBQM
from dimod.bqm.adjarraybqm import AdjArrayBQM

cdef cyAdjArrayBQM aabqm = AdjArrayBQM({'a': -1}, {'ab': 1}, 1.5, 'SPIN')

cdef cppAdjMapBQM[int, float] bqm = cppAdjMapBQM[int, float](aabqm.bqm_)


print(bqm.num_variables())
print(bqm.num_interactions())
