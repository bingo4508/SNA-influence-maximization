__author__ = 'bingo4508'

from networkx import nx
from model.LT import LTModel
from Strategy.util import *
import sys
import heapq as hq

run1_iter = 0
pr = None
NUM_PLAYERS = 2


def method(m, G, num_of_nodes):
    # Remove newly_activated/activated nodes
    for n in set(m.active_nodes_list):
        G.remove_node(n)
    for n in m.newly_active_nodes_list:
        G.remove_node(n)

    bc = nx.closeness_centrality(G, distance='influence')
    c = [(k, v) for k, v in bc.items()]
    c.sort(key=lambda tup: tup[1], reverse=True)

    selected_nodes = [c[i][0] for i in range(num_of_nodes)]
    return selected_nodes


def run():
    m = LTModel(nodes_file, edges_file, status_file, NUM_PLAYERS, player_id-1)

    g = m.G.reverse(copy=True)
    ns = method(m, g, nodes_num_per_iter)

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