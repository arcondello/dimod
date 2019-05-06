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
import unittest

import numpy as np

import dimod

from dimod.bqm.coo_bqm import CooBQM
from dimod.exceptions import WriteableError


class TestAggregate(unittest.TestCase):
    def test_duplicate(self):
        ldata = [0, 1]
        irow = [0, 0]
        icol = [1, 1]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        new = bqm.aggregate()

        np.testing.assert_array_equal(new.ldata, [0, 1])
        np.testing.assert_array_equal(new.irow, [0])
        np.testing.assert_array_equal(new.icol, [1])
        np.testing.assert_array_equal(new.qdata, [3])


class TestCooQuadratic(unittest.TestCase):
    def test_get_sorted_aggregated(self):
        ldata = [0, 1, 2]
        irow = [0, 1]
        icol = [1, 2]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        self.assertTrue(bqm.is_sorted)  # sanity check

        self.assertEqual(bqm.quadratic[(0, 1)], 1)
        self.assertEqual(bqm.quadratic[(1, 0)], 1)
        self.assertEqual(bqm.quadratic[(2, 1)], 2)
        self.assertEqual(bqm.quadratic[(1, 2)], 2)

    def test_get_unsorted_aggregated(self):
        ldata = [0, 1, 2]
        irow = [1, 1]
        icol = [0, 2]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        self.assertFalse(bqm.is_sorted)  # sanity check

        self.assertEqual(bqm.quadratic[(0, 1)], 1)
        self.assertEqual(bqm.quadratic[(1, 0)], 1)
        self.assertEqual(bqm.quadratic[(2, 1)], 2)
        self.assertEqual(bqm.quadratic[(1, 2)], 2)

    def test_get_sorted_duplicates(self):
        ldata = [0, 1, 2]
        irow = [0, 1, 1]
        icol = [1, 2, 2]
        qdata = [1, 2, 1]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        self.assertTrue(bqm.is_sorted)  # sanity check

        self.assertEqual(bqm.quadratic[(0, 1)], 1)
        self.assertEqual(bqm.quadratic[(1, 0)], 1)
        self.assertEqual(bqm.quadratic[(2, 1)], 3)
        self.assertEqual(bqm.quadratic[(1, 2)], 3)

    def test_get_unsorted_duplicates(self):
        ldata = [0, 1, 2]
        irow = [1, 1, 1]
        icol = [0, 2, 0]
        qdata = [1, 2, 3]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        self.assertFalse(bqm.is_sorted)  # sanity check

        self.assertEqual(bqm.quadratic[(0, 1)], 4)
        self.assertEqual(bqm.quadratic[(1, 0)], 4)
        self.assertEqual(bqm.quadratic[(2, 1)], 2)
        self.assertEqual(bqm.quadratic[(1, 2)], 2)

    def test_set_sorted_aggregated(self):
        ldata = [0, 1, 2]
        irow = [0, 1]
        icol = [1, 2]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        self.assertTrue(bqm.is_sorted)  # sanity check

        bqm.quadratic[(0, 1)] = 17

        np.testing.assert_array_equal(bqm.qdata, [17, 2])

        bqm.quadratic[(1, 0)] += 1

        np.testing.assert_array_equal(bqm.qdata, [18, 2])

        with self.assertRaises(KeyError):
            bqm.quadratic[('a', 0)] = 5

    # def test_get_unsorted_aggregated(self):
    #     ldata = [0, 1, 2]
    #     irow = [1, 1]
    #     icol = [0, 2]
    #     qdata = [1, 2]

    #     bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

    #     self.assertFalse(bqm.is_sorted)  # sanity check

    #     self.assertEqual(bqm.quadratic[(0, 1)], 1)
    #     self.assertEqual(bqm.quadratic[(1, 0)], 1)
    #     self.assertEqual(bqm.quadratic[(2, 1)], 2)
    #     self.assertEqual(bqm.quadratic[(1, 2)], 2)

    # def test_get_sorted_duplicates(self):
    #     ldata = [0, 1, 2]
    #     irow = [0, 1, 1]
    #     icol = [1, 2, 2]
    #     qdata = [1, 2, 1]

    #     bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

    #     self.assertTrue(bqm.is_sorted)  # sanity check

    #     self.assertEqual(bqm.quadratic[(0, 1)], 1)
    #     self.assertEqual(bqm.quadratic[(1, 0)], 1)
    #     self.assertEqual(bqm.quadratic[(2, 1)], 3)
    #     self.assertEqual(bqm.quadratic[(1, 2)], 3)

    # def test_get_unsorted_duplicates(self):
    #     ldata = [0, 1, 2]
    #     irow = [1, 1, 1]
    #     icol = [0, 2, 0]
    #     qdata = [1, 2, 3]

    #     bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

    #     self.assertFalse(bqm.is_sorted)  # sanity check

        # self.assertEqual(bqm.quadratic[(0, 1)], 4)
        # self.assertEqual(bqm.quadratic[(1, 0)], 4)
        # self.assertEqual(bqm.quadratic[(2, 1)], 2)
        # self.assertEqual(bqm.quadratic[(1, 2)], 2)


class TestDenseLinear(unittest.TestCase):
    def test_set_linear(self):
        bqm = CooBQM.from_ising({'b': 1}, {('a', 'b'): -2})

        bqm.linear['a'] = 5

        self.assertEqual(bqm.linear, {'a': 5, 'b': 1})


class TestFromIsing(unittest.TestCase):
    def test_empty(self):
        bqm = CooBQM.from_ising({}, {})

    def test_integer_labelled(self):
        bqm = CooBQM.from_ising({0: -1}, {(0, 1): 1})

        self.assertEqual(bqm.linear, {0: -1, 1: 0})
        self.assertEqual(bqm.quadratic, {(0, 1): 1})

    def test_str_labelled(self):
        bqm = CooBQM.from_ising({'a': -1}, {('a', 'b'): 1})

        self.assertEqual(bqm.linear, {'a': -1, 'b': 0})
        self.assertEqual(bqm.quadratic, {('a', 'b'): 1})

    def test_self_loop(self):
        bqm = CooBQM.from_ising({}, {('a', 'a'): -1})

        self.assertEqual(bqm.linear, {'a': 0})
        self.assertEqual(bqm.quadratic, {('a', 'a'): -1})  # no reduction


class TestSorted(unittest.TestCase):
    def test_duplicate(self):
        # nothing should change
        ldata = [0, 1]
        irow = [0, 0]
        icol = [1, 1]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        self.assertTrue(bqm.is_sorted)

        bqm.sort()

        np.testing.assert_array_equal(bqm.irow, irow)
        np.testing.assert_array_equal(bqm.icol, icol)
        np.testing.assert_array_equal(bqm.qdata, qdata)

    def test_unsorted(self):
        ldata = [0, 1, 2]
        irow = [0, 0]
        icol = [2, 1]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        self.assertFalse(bqm.is_sorted)

        bqm.sort()

        np.testing.assert_array_equal(bqm.irow, [0, 0])
        np.testing.assert_array_equal(bqm.icol, [1, 2])
        np.testing.assert_array_equal(bqm.qdata, [2, 1])

    def test_immutable_sorted(self):
        ldata = [0, 1]
        irow = [0, 0]
        icol = [1, 1]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        bqm.setflags(write=False)

        bqm.sort()  # already sorted so no error

    def test_immutable_unsorted(self):
        ldata = [0, 1, 2]
        irow = [0, 0]
        icol = [2, 1]
        qdata = [1, 2]

        bqm = CooBQM(ldata, (irow, icol, qdata), 0, 'SPIN')

        bqm.setflags(write=False)

        with self.assertRaises(WriteableError) as e:
            bqm.sort()

        # make sure the type and attr are correct in the error message
        msg = e.exception.args[0]
        self.assertEqual(msg,
                         ("cannot be sorted while {}.{} is set to "
                          'False'.format(type(bqm).__name__, 'is_writeable')))


class TestWriteable(unittest.TestCase):
    def test_setflags(self):
        bqm = CooBQM.from_ising({'a': -1}, {('a', 'b'): 1})
        self.assertTrue(bqm.is_writeable)

        bqm.relabel({'a': 0})

        self.assertEqual(bqm.linear, {0: -1, 'b': 0})
        self.assertEqual(bqm.quadratic, {(0, 'b'): 1})

        bqm.setflags(write=False)
        self.assertFalse(bqm.is_writeable)

        with self.assertRaises(ValueError):
            bqm.relabel({'b': 1})

        bqm.setflags(write=True)
        self.assertTrue(bqm.is_writeable)

        bqm.relabel({'b': 1})

        self.assertEqual(bqm.linear, {0: -1, 1: 0})
        self.assertEqual(bqm.quadratic, {(0, 1): 1})

    def test_set_linear(self):
        bqm = CooBQM.from_ising({'b': 1}, {('a', 'b'): -2})
        bqm.setflags(write=False)

        with self.assertRaises(WriteableError) as e:
            bqm.linear['a'] = 5

        # make sure the type and attr are correct in the error message
        msg = e.exception.args[0]
        self.assertEqual(msg,
                         ('Cannot set linear bias when {}.{} is '
                          'False'.format(type(bqm).__name__, 'is_writeable')))