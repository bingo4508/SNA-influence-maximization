__author__ = 'bingo4508'

from networkx import nx
from model.LT import LTModel
from Strategy.util import *
import sys

run1_iter = 0
pr = None
NUM_PLAYERS = 2

# Graph G need to be reversed in advance
# Run page rank only once
def method1(G, num_of_nodes, alpha=0.85, weight=False):
    global run1_iter, pr
    if run1_iter == 0:
        if weight:
            pr = nx.pagerank(G, weight='influence', alpha=alpha)
        else:
            pr = nx.pagerank(G, alpha=alpha)
    node_rank = [(k, v) for k, v in pr.items()]
    node_rank.sort(key=lambda tup: tup[1], reverse=True)

    nodes = [node_rank[run1_iter*num_of_nodes+j][0] for j in range(num_of_nodes)]
    run1_iter += 1
    return nodes


# Run page rank every iteration
def method2(model, G, num_of_nodes, alpha=0.85, weight=False):
    # Remove newly_activated/activated nodes
    for n in set(model.active_nodes_list):
        G.remove_node(n)
    for n in model.newly_active_nodes_list:
        G.remove_node(n)

    if weight:
        pr = nx.pagerank(G, weight='influence', alpha=alpha)
    else:
        pr = nx.pagerank(G, alpha=alpha)

    node_rank = [(k, v) for k, v in pr.items()]
    node_rank.sort(key=lambda tup: tup[1], reverse=True)
    nodes = [node_rank[i][0] for i in range(num_of_nodes)]
    return nodes


def run():
    m = LTModel(nodes_file, edges_file, status_file, NUM_PLAYERS, player_id-1)
    # Page Rank - run once, performance bad
    # ns = run1(g, PICK_NODE, 0.9, weight=True)

    # Page Rank - run every iteration
    g = m.G.reverse(copy=True)
    ns = method2(m, g, nodes_num_per_iter, 0.9, weight=True)

    # Write output
    write_nodes_list(ns, selected_nodes_file)


def test(p_player_id, p_nodes_file, p_edges_file, p_status_file, p_selected_nodes_file, p_nodes_num_per_iter, p_time_limit_in_sec):
    global player_id, nodes_file, edges_file, status_file, nodes_num_per_iter, selected_nodes_file, time_limit_in_sec

    player_id = int(p_player_id)
    nodes_file = p_nodes_file
    edges_file = p_edges_file
    status_file = p_status_file
    nodes_num_per_iter = int(p_nodes_num_per_iter)
    selected_nodes_file = p_selected_nodes_file
    time_limit_in_sec = float(p_time_limit_in_sec)

    run()

# main function
# if __name__ == '__main__':
#     player_id = int(sys.argv[1])
#     nodes_file = sys.argv[2]
#     edges_file = sys.argv[3]
#     status_file = sys.argv[4]
#     nodes_num_per_iter = int(sys.argv[5])
#     selected_nodes_file = sys.argv[6]
#     time_limit_in_sec = int(sys.argv[7])
#     run()