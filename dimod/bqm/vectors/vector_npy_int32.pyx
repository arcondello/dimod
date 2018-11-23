# distutils: language = c++
# cython: language_level = 3
#
# NOTE: This is a procedurally generated file. It should not be edited. See vector.pyx.template
#

import collections.abc as abc

from libcpp.vector cimport vector

import numpy as np
cimport numpy as np


ctypedef np.npy_int32 dtype


cdef class _Vector_npy_int32:
    cdef vector[dtype] biases

    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]

    def __cinit__(self):
        self.biases = vector[dtype]()

    def __init__(self, iterable=tuple()):
        for v in iterable:
            self.biases.push_back(v)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t itemsize = sizeof(self.biases[0])

        self.shape[0] = self.biases.size()

        self.strides[0] = sizeof(self.biases[0])

        buffer.buf = <char *>&(self.biases[0])
        buffer.format = 'i'  # needs to correspond to dtype
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = self.biases.size() * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    def __len__(self):
        return self.biases.size()

    def __getitem__(self, Py_ssize_t i):
        if i < 0 or i >= self.biases.size():
            raise IndexError('index out of range')
        return self.biases[i]

    def __delitem__(self, Py_ssize_t i):
        if i < 0 or i >= self.biases.size():
            raise IndexError('assignment index out of range')
        self.biases.erase(self.biases.begin() + i)

    def __setitem__(self, Py_ssize_t i, dtype bias):
        if i < 0 or i >= self.biases.size():
            raise IndexError('assignment index out of range')
        self.biases[i] = bias

    def insert(self, Py_ssize_t i, dtype bias):
        if i >= len(self):
            self.biases.push_back(bias)
        else:
            self.biases[i] = bias

class Vector_npy_int32(_Vector_npy_int32, abc.MutableSequence):
    __slots__ = ()

    def __str__(self):
        return str(list(self))

    def __repr__(self):
        return 'vector(%s, dtype=%r' % (self, np.asarray(self).dtype.name)
