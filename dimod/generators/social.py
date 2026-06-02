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

import dimod

__all__ = ["structural_imbalance",
           ]


def structural_imbalance(graph: 'nx.Graph') -> dimod.BinaryQuadraticModel:
    """Construct a binary quadratic model (BQM) to calculate the structural
    imbalance of a signed social network.

    A signed social network graph is a graph whose signed edges represent
    friendly/hostile interactions between nodes. A signed social network is
    considered balanced if it can be cleanly divided into two factions, where
    all relations within a faction are friendly, and all relations between
    factions are hostile. The measure of imbalance or frustration is the minimum
    number of edges that violate this rule.

    Args:
        graph:
            A social graph (in a NetworkX graph) on which each edge has a 'sign'
            attribute with a numeric value.

    Returns:
        A binary quadratic model. Each variable in the model represents a node
        in the signed social network. The solution that minimized the BQM will
        assign each variable a value, either -1 or 1. This bi-coloring defines
        the factions.

    Raises:
        ValueError: If any edge does not have a 'sign' attribute.

    Examples:
        >>> import dimod
        >>> import networkx as nx
        ...
        >>> S = nx.Graph()
        >>> S.add_edge('Alice', 'Bob', sign=1)  # Alice and Bob are friendly
        >>> S.add_edge('Alice', 'Eve', sign=-1)  # Alice and Eve are hostile
        >>> S.add_edge('Bob', 'Eve', sign=-1)  # Bob and Eve are hostile
        ...
        >>> bqm = dimod.generators.structural_imbalance(S)
        >>> bqm.linear      # doctest: +SKIP
        {'Alice': 0.0, 'Bob': 0.0, 'Eve': 0.0}
        >>> bqm.quadratic   # doctest: +SKIP
        {('Alice', 'Bob'): -1.0, ('Alice', 'Eve'): 1.0, ('Bob', 'Eve'): 1.0}

    """
    if not hasattr(graph, 'nodes') or not hasattr(graph, 'edges'):
        raise ValueError("Signed social network graph in NetworkX format required")

    bqm = dimod.BinaryQuadraticModel.empty('SPIN')

    for v in graph:
        bqm.add_linear(v, 0.0)

    for u, v, data in graph.edges(data=True):
        try:
            bqm.add_quadratic(u, v, -1. * data['sign'])
        except KeyError:
            raise ValueError("graph should be a signed social graph, "
                             "each edge should have a 'sign' attr")

    return bqm
