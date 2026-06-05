# Copyright 2026 D-Wave
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

import itertools
import math

import dimod
from dimod.typing import GraphLike
from dimod.decorators import graph_argument

__all__ = ["graph_partition",
           ]


@graph_argument('graph', as_networkx=True)
def graph_partition(graph: GraphLike,
                    num_partitions: int,
                    ) -> dimod.ConstrainedQuadraticModel:
    """Find a constrained quadratic model for the graph's partitions.

    Defines a CQM with ground states corresponding to a balanced k-partition
    of ``graph``. A k-partition is a collection of k subsets of the vertices
    of ``graph`` such that each vertex is in exactly one subset, and the number
    of edges between vertices in different subsets is as small as possible.
    If ``graph`` is a weighted graph, the sum of weights over those edges are
    minimized.

    Args:
        graph:
            The graph to partition. Either an integer ``n``, interpreted as a
            complete graph of size ``n``, a nodes/edges pair, a list of edges or
            a NetworkX graph. When NetworkX graph is provided, optional edge
            weights can be provided in the ``weight`` attribute.

        num_partitions:
            The number of subsets in the desired partition.

    Returns:
        A constrained quadratic model with ground states corresponding to a
        partition problem. The nodes of ``graph`` are discrete logical variables
        of the CQM, where the cases are the different partitions the node
        can be assigned to. The objective is given as the number of edges
        connecting nodes in different partitions.

    """
    partition_size = graph.number_of_nodes() / num_partitions
    partitions = range(num_partitions)
    cqm = dimod.ConstrainedQuadraticModel()

    # Variables will be added using the discrete method in CQM
    x = {vk: dimod.Binary(vk) for vk in itertools.product(graph.nodes, partitions)}

    for v in graph.nodes:
        cqm.add_discrete(((v, k) for k in partitions), label=v)

    if not math.isclose(partition_size, int(partition_size)):
        # if number of nodes don't divide into num_partitions,
        # accept partitions of size ceil() or floor()
        floor, ceil = int(partition_size), int(partition_size+1)
        for k in partitions:
            cqm.add_constraint(
                dimod.quicksum((x[u, k] for u in graph.nodes)) >= floor,
                label=f'equal_partition_low_{k}')
            cqm.add_constraint(
                dimod.quicksum((x[u, k] for u in graph.nodes)) <= ceil,
                label=f'equal_partition_high_{k}')
    else:
        # each partition must have partition_size elements
        for k in partitions:
            cqm.add_constraint(
                dimod.quicksum((x[u, k] for u in graph.nodes)) == int(partition_size),
                label=f'equal_partition_{k}')

    cuts = 0
    for (u, v, d) in graph.edges(data=True):
        for k in partitions:
            w = d.get('weight',1)
            cuts += w * x[u,k] * x[v,k]

    if cuts:
        cqm.set_objective(-cuts)

    return cqm
