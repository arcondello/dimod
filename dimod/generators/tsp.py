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

import dimod

__all__ = ["traveling_salesperson",
           "traveling_salesman",
           ]


def traveling_salesperson(graph: 'nx.Graph',
                          lagrange: float | None = None,
                          weight: str = 'weight',
                          missing_edge_weight: float | None = None
                          ) -> dimod.BinaryQuadraticModel:
    """Return a binary quadratic model (BQM) with ground states corresponding
    to a minimum TSP route.

    If :math:`|graph|` is the number of nodes in the graph, the resulting qubo
    will have:

    * :math:`|graph|^2` variables/nodes
    * :math:`2 |graph|^2 (|graph| - 1)` interactions/edges

    Args:
        graph:
            A complete graph in which each edge has an attribute giving its
            weight, given as a NetworkX graph.

        lagrange:
            Lagrange parameter to weight constraints (no edges within set)
            versus objective (largest set possible).

        weight:
            The name of the edge attribute containing the weight, defaults to
            "weight".
    
        missing_edge_weight:
            For bi-directional graphs, the weight given to missing edges.
            If None is given (the default), missing edges will be set to
            the sum of all weights.

    Returns:
       A binary quadratic model (BQM) with ground states corresponding to a
       minimum traveling salesperson route. The BQM variables are labelled
       ``(c, t)`` where ``c`` is a node in ``graph`` and ``t`` is the time
       index. For instance, if ``('a', 0)`` is 1 in the ground state, that means
       the node 'a' is visited first.

    """
    if not hasattr(graph, 'nodes') or not hasattr(graph, 'edges'):
        raise ValueError("A NetworkX graph with weights data required")

    N = graph.number_of_nodes()

    if lagrange is None:
        # If no lagrange parameter provided, set to 'average' tour length.
        # Usually a good estimate for a lagrange parameter is between 75-150%
        # of the objective function value, so we come up with an estimate for 
        # tour length and use that.
        if graph.number_of_edges() > 0:
            lagrange = (graph.size(weight=weight)
                        * graph.number_of_nodes()
                        / graph.number_of_edges())
        else:
            lagrange = 2

    # calculate default missing_edge_weight if required
    if missing_edge_weight is None:
        # networkx method to calculate sum of all weights
        missing_edge_weight = graph.size(weight=weight)

    # some input checking
    if N in (1, 2):
        msg = "graph must have at least 3 nodes or be empty"
        raise ValueError(msg)

    # Creating the BQM
    bqm = dimod.BinaryQuadraticModel.empty(dimod.BINARY)

    # Constraint that each row has exactly one 1
    for node in graph:
        for pos_1 in range(N):
            bqm.add_linear((node, pos_1), -lagrange)
            for pos_2 in range(pos_1+1, N):
                bqm.add_quadratic((node, pos_1), (node, pos_2), 2.0*lagrange)

    # Constraint that each col has exactly one 1
    for pos in range(N):
        for node_1 in graph:
            bqm.add_linear((node_1, pos), -lagrange)
            for node_2 in set(graph)-{node_1}:
                # quadratic coefficient is 2*lagrange, but we are placing this value 
                # above *and* below the diagonal, so we put half in each position.
                bqm.add_quadratic((node_1, pos), (node_2, pos), lagrange)

    # Objective that minimizes distance
    for u, v in itertools.combinations(graph.nodes, 2):
        for pos in range(N):
            nextpos = (pos + 1) % N

            # going from u -> v
            try:
                value = graph[u][v][weight]
            except KeyError:
                value = missing_edge_weight

            bqm.add_quadratic((u, pos), (v, nextpos), value)

            # going from v -> u
            try:
                value = graph[v][u][weight]
            except KeyError:
                value = missing_edge_weight

            bqm.add_quadratic((v, pos), (u, nextpos), value)

    return bqm


traveling_salesman = traveling_salesperson
