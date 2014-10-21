__author__ = 'bingo4508'

from networkx import nx
from model.LT import LTModel
from Strategy.util import *
import sys

run1_iter = 0
pr = None
NUM_PLAYERS = 2


def method(graph, num_of_nodes):
    # write your code here
    nodes_list = graph.nodes(data=True)
    candidate_list = list()

    for n in nodes_list:
        if n[1]['status'] == 'inactivated':
            candidate_list.append((n[0], graph.out_degree(n[0])))

    candidate_list.sort(key = lambda x: x[1], reverse = True)

    # store the selected nodes into return_nodes_list
    return_num = min(num_of_nodes, len(candidate_list))
    return_nodes_list = list()
    for i in range(0, return_num):
        return_nodes_list.append(candidate_list[i][0])

    return return_nodes_list


def run():
    m = LTModel(nodes_file, edges_file, status_file, NUM_PLAYERS, player_id-1)

    g = m.get_copy_graph()
    ns = method(g, nodes_num_per_iter)

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