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

import unittest

import dimod

from dimod import BQM, Spin, Binary, CQM, Integer
from dimod.sym import Sense


class TestAddConstraint(unittest.TestCase):
    def test_bqm(self):
        cqm = CQM()

        bqm = BQM({'a': -1}, {'ab': 1}, 1.5, 'SPIN')
        cqm.add_constraint(bqm, '<=')
        cqm.add_constraint(bqm, '>=')  # add it again

    def test_duplicate(self):
        cqm = CQM()

        bqm = BQM({'a': -1}, {'ab': 1}, 1.5, 'SPIN')

        cqm.add_constraint(bqm <= 5, label='hello')
        with self.assertRaises(ValueError):
            cqm.add_constraint(bqm <= 5, label='hello')

    def test_symbolic(self):
        cqm = CQM()

        a = Spin('a')
        b = Spin('b')
        c = Spin('c')

        cqm.add_constraint(2*a*b + b*c - c + 1 <= 1)

    def test_symbolic_mixed(self):
        cqm = CQM()

        x = Binary('x')
        s = Spin('s')
        i = Integer('i')

        cqm.add_constraint(2*i + s + x <= 2)

        self.assertIs(cqm.vartype('x'), dimod.BINARY)
        self.assertIs(cqm.vartype('s'), dimod.SPIN)
        self.assertIs(cqm.vartype('i'), dimod.INTEGER)

    def test_terms(self):
        cqm = CQM()

        a = cqm.add_variable('a', 'BINARY')
        b = cqm.add_variable('b', 'BINARY')
        c = cqm.add_variable('c', 'INTEGER')

        cqm.add_constraint([(a, b, 1), (b, 2.5,), (3,), (c, 1.5)], sense='<=')

    def test_terms_integer_bounds(self):
        # bug report: https://github.com/dwavesystems/dimod/issues/943
        cqm = CQM()
        i = Integer('i', lower_bound=-1, upper_bound=5)
        cqm.set_objective(i)
        label = cqm.add_constraint([('i', 1)], sense='<=')  # failing in #943

        self.assertEqual(cqm.constraints[label].lhs.lower_bound('i'), -1)
        self.assertEqual(cqm.constraints[label].lhs.upper_bound('i'), 5)


class TestAddDiscrete(unittest.TestCase):
    def test_simple(self):
        cqm = CQM()
        cqm.add_discrete('abc')

    def test_label(self):
        cqm = CQM()
        label = cqm.add_discrete('abc', label='hello')
        self.assertEqual(label, 'hello')
        self.assertEqual(cqm.variables, 'abc')


class TestAdjVector(unittest.TestCase):
    # this will be deprecated in the future
    def test_construction(self):
        cqm = CQM()

        with self.assertWarns(DeprecationWarning):
            bqm = dimod.AdjVectorBQM({'ab': 1}, 'SPIN')

        cqm.set_objective(bqm)
        label = cqm.add_constraint(bqm, sense='==', rhs=1)
        self.assertIsInstance(cqm.objective, dimod.QuadraticModel)
        self.assertIsInstance(cqm.constraints[label].lhs, BQM)


class TestBounds(unittest.TestCase):
    def test_inconsistent(self):
        i0 = Integer('i')
        i1 = Integer('i', upper_bound=1)
        i2 = Integer('i', lower_bound=-2)

        cqm = CQM()
        cqm.set_objective(i0)
        with self.assertRaises(ValueError):
            cqm.add_constraint(i1 <= 1)

        cqm = CQM()
        cqm.add_constraint(i0 <= 1)
        with self.assertRaises(ValueError):
            cqm.set_objective(i2)

    def test_later_defn(self):
        i0 = Integer('i')
        i1 = Integer('i', upper_bound=1)

        cqm = CQM()
        cqm.add_variable('i', 'INTEGER')
        cqm.set_objective(i0)
        with self.assertRaises(ValueError):
            cqm.add_constraint(i1 <= 1)

        cqm.add_variable('i', 'INTEGER')


class TestCopy(unittest.TestCase):
    def test_deepcopy(self):
        from copy import deepcopy

        cqm = CQM()

        x = Binary('x')
        s = Spin('s')
        i = Integer('i')

        cqm.set_objective(i + s + x)
        constraint = cqm.add_constraint(i + s + x <= 1)

        new = deepcopy(cqm)

        self.assertTrue(new.objective.is_equal(cqm.objective))
        self.assertTrue(new.constraints[constraint].lhs.is_equal(cqm.constraints[constraint].lhs))


class TestCQMtoBQM(unittest.TestCase):
    def test_empty(self):
        bqm, inverter = dimod.cqm_to_bqm(dimod.CQM())
        self.assertEqual(bqm.shape, (0, 0))
        self.assertEqual(bqm.vartype, dimod.BINARY)

    def test_bqm_objective_only(self):
        x, y, z = dimod.Binaries('xyz')

        cqm = CQM.from_bqm(x*y + 2*y*z + 8*x + 5)

        bqm, inverter = dimod.cqm_to_bqm(cqm)

        self.assertEqual(bqm, x*y + 2*y*z + 8*x + 5)

    def test_qm_objective_only(self):
        i = dimod.Integer('i', upper_bound=7)
        j = dimod.Integer('j', upper_bound=9)
        x = dimod.Binary('x')

        qm = i*j + 5*j*x + 8*i + 3*x + 5
        cqm = CQM.from_qm(qm)

        bqm, inverter = dimod.cqm_to_bqm(cqm)

        sampleset = dimod.ExactSolver().sample(bqm)

        for bin_sample, energy in sampleset.data(['sample', 'energy']):
            int_sample = inverter(bin_sample)
            self.assertEqual(qm.energy(int_sample), energy)

    def test_bqm_equality_constraint(self):
        x, y, z = dimod.Binaries('xyz')

        cqm = CQM()
        cqm.add_constraint(x + y + z == 1)

        bqm, inverter = dimod.cqm_to_bqm(cqm, lagrange_multiplier=1)

        self.assertEqual(bqm, (x + y + z - 1)*(x + y + z - 1))

    def test_bqm_equality_constraint_offset(self):
        x, y, z = dimod.Binaries('xyz')

        cqm = CQM()
        cqm.add_constraint(x + y + z - 1 == 2)

        bqm, inverter = dimod.cqm_to_bqm(cqm, lagrange_multiplier=1)

        self.assertEqual(bqm, (x + y + z - 3)*(x + y + z - 3))

    def test_bqm_Le_constraint(self):
        x, y, z = dimod.Binaries('xyz')

        cqm = CQM()
        cqm.add_constraint(x + y + z <= 2)

        bqm, inverter = dimod.cqm_to_bqm(cqm, lagrange_multiplier=1)

        configs = set()
        for sample in dimod.ExactSolver().sample(bqm).lowest().samples():
            self.assertLessEqual(sample['x'] + sample['y'] + sample['z'], 2)
            configs.add((sample['x'], sample['y'], sample['z']))
        self.assertEqual(len(configs), 7)

    def test_bqm_Ge_constraint(self):
        x, y, z = dimod.Binaries('xyz')

        cqm = CQM()
        cqm.add_constraint(x + y + z >= 2)

        bqm, inverter = dimod.cqm_to_bqm(cqm, lagrange_multiplier=1)

        configs = set()
        for sample in dimod.ExactSolver().sample(bqm).lowest().samples():
            self.assertGreaterEqual(sample['x'] + sample['y'] + sample['z'], 2)
            configs.add((sample['x'], sample['y'], sample['z']))
        self.assertEqual(len(configs), 4)

    def test_qm_Ge_constraint(self):
        i = dimod.Integer('i', upper_bound=7)
        j = dimod.Integer('j', upper_bound=9)
        x = dimod.Binary('x')

        cqm = CQM()
        cqm.add_constraint(i + j + x >= 5)

        bqm, inverter = dimod.cqm_to_bqm(cqm, lagrange_multiplier=1)

        for bin_sample in dimod.ExactSolver().sample(bqm).lowest().samples():
            int_sample = inverter(bin_sample)
            self.assertGreaterEqual(int_sample['i'] + int_sample['j'] + int_sample['x'], 5)


class TestDeprecated(unittest.TestCase):
    def test_typedvariables(self):
        cqm = CQM()

        x = Binary('x')
        s = Spin('s')
        i = Integer('i')

        cqm.set_objective(x + i)
        cqm.add_constraint(i + s <= 1)

        with self.assertWarns(DeprecationWarning):
            self.assertEqual(cqm.variables.lower_bounds, {'x': 0.0, 'i': 0.0, 's': -1.0})
        with self.assertWarns(DeprecationWarning):
            self.assertEqual(cqm.variables.upper_bounds, {'x': 1.0, 'i': 9007199254740991.0, 's': 1})
        with self.assertWarns(DeprecationWarning):
            self.assertEqual(len(cqm.variables.lower_bounds), 3)
        with self.assertWarns(DeprecationWarning):
            self.assertEqual(len(cqm.variables.upper_bounds), 3)

        with self.assertWarns(DeprecationWarning):
            self.assertEqual(list(cqm.variables.vartypes),
                             [dimod.BINARY, dimod.INTEGER, dimod.SPIN])

        with self.assertWarns(DeprecationWarning):
            self.assertIs(cqm.variables.vartype('x'), dimod.BINARY)


class TestFromDQM(unittest.TestCase):
    def test_case_label(self):
        dqm = dimod.CaseLabelDQM()
        u = dqm.add_variable({'red', 'green', 'blue'}, shared_labels=True)
        v = dqm.add_variable(['blue', 'yellow', 'brown'], label='v', shared_labels=True)
        dqm.set_linear_case(u, 'red', 1)
        dqm.set_linear_case(v, 'yellow', 2)
        dqm.set_quadratic_case(u, 'green', v, 'blue', -0.5)
        dqm.set_quadratic_case(u, 'blue', v, 'brown', -0.5)

        cqm = CQM.from_discrete_quadratic_model(dqm)

        self.assertEqual(cqm.objective.linear,
                         {(0, 'blue'): 0.0, (0, 'red'): 1.0, (0, 'green'): 0.0,
                          ('v', 'blue'): 0.0, ('v', 'brown'): 0.0, ('v', 'yellow'): 2.0})
        self.assertEqual(cqm.objective.quadratic,
                         {(('v', 'blue'), (0, 'green')): -0.5, (('v', 'brown'), (0, 'blue')): -0.5})
        self.assertEqual(cqm.objective.offset, 0)
        self.assertTrue(all(cqm.objective.vartype(v) is dimod.BINARY
                            for v in cqm.objective.variables))

        self.assertEqual(set(cqm.constraints), set(dqm.variables))

    def test_empty(self):
        dqm = dimod.DiscreteQuadraticModel()

        cqm = CQM.from_discrete_quadratic_model(dqm)

        self.assertEqual(len(cqm.variables), 0)
        self.assertEqual(len(cqm.constraints), 0)

    def test_typical(self):
        dqm = dimod.DQM()
        u = dqm.add_variable(4)
        v = dqm.add_variable(3)
        dqm.set_quadratic(u, v, {(0, 2): -1, (2, 1): 1})
        dqm.set_linear(u, [0, 1, 2, 3])
        dqm.offset = 5

        cqm = CQM.from_discrete_quadratic_model(dqm)

        self.assertEqual(cqm.variables, [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)])

        self.assertEqual(cqm.objective.linear,
                         {(0, 0): 0, (0, 1): 1, (0, 2): 2, (0, 3): 3,
                          (1, 1): 0, (1, 2): 0, (1, 0): 0})
        self.assertEqual(cqm.objective.quadratic, {((1, 1), (0, 2)): 1.0, ((1, 2), (0, 0)): -1.0})
        self.assertEqual(cqm.objective.offset, dqm.offset)

        # keys of constraints are the variables of DQM
        self.assertEqual(set(cqm.constraints), set(dqm.variables))


class FromQM(unittest.TestCase):
    def test_from_bqm(self):
        bqm = BQM({'a': -1}, {'ab': 1}, 1.5, 'SPIN')
        cqm = CQM.from_bqm(bqm)
        self.assertTrue(cqm.objective.is_equal(dimod.QM.from_bqm(bqm)))

    def test_from_qm(self):
        qm = BQM({'a': -1}, {'ab': 1}, 1.5, 'SPIN') + Integer('i')
        cqm = CQM.from_quadratic_model(qm)
        self.assertTrue(cqm.objective.is_equal(qm))


class TestIterViolations(unittest.TestCase):
    def test_no_constraints(self):
        cqm = CQM.from_bqm(-Binary('a') + Binary('a')*Binary('b') + 1.5)
        self.assertEqual(cqm.violations({'a': 1, 'b': 1}), {})

    def test_binary(self):
        a, b, c = dimod.Binaries('abc')

        cqm = CQM()
        cqm.add_constraint(a + b + c == 1, label='onehot')
        cqm.add_constraint(a*b <= 0, label='ab LE')
        cqm.add_constraint(c >= 1, label='c GE')

        sample = {'a': 0, 'b': 0, 'c': 1}  # satisfying sample
        self.assertEqual(cqm.violations(sample), {'ab LE': 0.0, 'c GE': 0.0, 'onehot': 0.0})
        self.assertEqual(cqm.violations(sample, skip_satisfied=True), {})

        sample = {'a': 1, 'b': 0, 'c': 0}  # violates one
        self.assertEqual(cqm.violations(sample), {'ab LE': 0.0, 'c GE': 1.0, 'onehot': 0.0})
        self.assertEqual(cqm.violations(sample, skip_satisfied=True), {'c GE': 1.0})

    def test_integer(self):
        i = dimod.Integer('i', lower_bound=-1000)
        j, k = dimod.Integers('jk')

        cqm = CQM()
        label_le = cqm.add_constraint(i + j*k <= 5)
        label_ge = cqm.add_constraint(i + j >= 1000)

        sample = {'i': 105, 'j': 4, 'k': 5}
        self.assertEqual(cqm.violations(sample), {label_le: 120.0, label_ge: 891.0})


class TestSerialization(unittest.TestCase):
    def test_functional(self):
        cqm = CQM()

        bqm = BQM({'a': -1}, {'ab': 1}, 1.5, 'SPIN')
        cqm.add_constraint(bqm, '<=')
        cqm.add_constraint(bqm, '>=')
        cqm.set_objective(BQM({'c': -1}, {}, 'SPIN'))
        cqm.add_constraint(Spin('a')*Integer('d')*5 <= 3)

        new = CQM.from_file(cqm.to_file())

        self.assertTrue(cqm.objective.is_equal(new.objective))
        self.assertEqual(set(cqm.constraints), set(new.constraints))
        for label, constraint in cqm.constraints.items():
            self.assertTrue(constraint.lhs.is_equal(new.constraints[label].lhs))
            self.assertEqual(constraint.rhs, new.constraints[label].rhs)
            self.assertEqual(constraint.sense, new.constraints[label].sense)

    def test_functional_empty(self):
        new = CQM.from_file(CQM().to_file())

    def test_functional_discrete(self):
        cqm = CQM()

        bqm = BQM({'a': -1}, {'ab': 1}, 1.5, 'SPIN')
        cqm.add_constraint(bqm, '<=')
        cqm.add_constraint(bqm, '>=')
        cqm.set_objective(Integer('c'))
        cqm.add_constraint(Spin('a')*Integer('d')*5 <= 3)
        cqm.add_discrete('efg')

        new = CQM.from_file(cqm.to_file())

        self.assertTrue(cqm.objective.is_equal(new.objective))
        self.assertEqual(set(cqm.constraints), set(new.constraints))
        for label, constraint in cqm.constraints.items():
            self.assertTrue(constraint.lhs.is_equal(new.constraints[label].lhs))
            self.assertEqual(constraint.rhs, new.constraints[label].rhs)
            self.assertEqual(constraint.sense, new.constraints[label].sense)
        self.assertSetEqual(cqm.discrete, new.discrete)

    def test_header(self):
        from dimod.serialization.fileview import read_header

        cqm = CQM()

        x = Binary('x')
        s = Spin('s')
        i = Integer('i')

        cqm.set_objective(x + 3*i + s*x)
        cqm.add_constraint(x*s + x <= 5)
        cqm.add_constraint(i*i + i*s <= 4)

        header_info = read_header(cqm.to_file(), b'DIMODCQM')

        self.assertEqual(header_info.data,
                         {'num_biases': 11, 'num_constraints': 2,
                          'num_quadratic_variables': 4, 'num_variables': 3})


class TestSetObjective(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(CQM().objective.num_variables, 0)

    def test_set(self):
        cqm = CQM()
        cqm.set_objective(Integer('a') * 5)
        self.assertTrue(cqm.objective.is_equal(Integer('a') * 5))


class TestSubstituteSelfLoops(unittest.TestCase):
    def test_typical(self):
        i, j, k = dimod.Integers('ijk')

        cqm = CQM()

        cqm.set_objective(2*i*i + i*k)
        label = cqm.add_constraint(-3*i*i + 4*j*j <= 5)

        mapping = cqm.substitute_self_loops()

        self.assertIn('i', mapping)
        self.assertIn('j', mapping)
        self.assertEqual(len(mapping), 2)

        self.assertEqual(cqm.objective.quadratic, {('k', 'i'): 1, (mapping['i'], 'i'): 2})
        self.assertEqual(cqm.constraints[label].lhs.quadratic,
                         {(mapping['i'], 'i'): -3.0, (mapping['j'], 'j'): 4.0})
        self.assertEqual(len(cqm.constraints), 3)

        for v, new in mapping.items():
            self.assertIn(new, cqm.constraints)
            self.assertEqual(cqm.constraints[new].sense, Sense.Eq)
            self.assertEqual(cqm.constraints[new].lhs.linear, {v: 1, new: -1})
            self.assertEqual(cqm.constraints[new].rhs, 0)
