from Strategy import page_rank as pr
import time

__author__ = 'bingo4508'

from model.LT import LTModel
import Strategy.page_rank as pr
import max_weight_strategy as mws
import max_outdegree_strategy as mos
import Strategy.lazy_greedy as lg
import Strategy.centrality as c

# m = LTModel('../networks/egofb_lt_nodes.txt', '../networks/egofb_lt_edges.txt', 'status.txt', 2, 1)
m = LTModel('../networks/hepth_lt_nodes.txt', '../networks/hepth_lt_edges.txt', 'status.txt', 2, 1)


PICK_NODE = 10
ITERATION = 10

for i in range(ITERATION):
    t_start = time.clock()

    ##################################################################
    # Page Rank - run once
    # ns = pr.method1(g, PICK_NODE, 0.9)
    # ns = pr.method1(g, PICK_NODE, 0.9, weight=True)

    # Page Rank - run every iteration
    # g = m.G.reverse(copy=True)
    # ns = pr.method2(m, g, PICK_NODE, 0.9)
    # ns = pr.method2(m, g, PICK_NODE, 0.9, weight=True)

    # Max Weight
    # ns = mws.method(m.G, PICK_NODE)

    # Max Out degree
    # ns = mos.run(m.G, PICK_NODE)

    # Lazy Greedy
    ns = lg.method(m, PICK_NODE)

    # Betweenness centrality
    # g = m.G.reverse(copy=True)
    # ns = c.method(m, g, PICK_NODE)

    ##################################################################
    print("Time: %f secs" % (time.clock()-t_start))
    # Propagation...
    m.select_nodes(ns)
    m.propagate()
    #
    num_active = m.get_activated_nodes()[m.my_id] - (PICK_NODE*(i+1))
    print("Activated: %d" % num_active)