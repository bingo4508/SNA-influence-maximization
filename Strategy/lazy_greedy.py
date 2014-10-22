__author__ = 'bingo4508'

from networkx import nx
from model.LT import LTModel
from Strategy.util import *
import sys
import heapq as hq

run1_iter = 0
pr = None
NUM_PLAYERS = 2


def calculate_spread(m, node_list):
    ori_active = m.get_activated_nodes()[m.my_id]
    m.select_nodes(node_list, test=True)
    m.propagate(test=True)
    spread = m.get_activated_nodes()[m.my_id]
    m.reset()
    return spread-ori_active


def method(m, num_of_nodes):
    selected_nodes = []
    q_table = []    # (margin, node, iter)
    iteration = 1
    for n in m.G.nodes():
        if m.G.node[n]['status'] == 'inactivated':
            margin = -1*calculate_spread(m, [n])
        else:
            margin = 0
        hq.heappush(q_table, (margin, n, 1))

    while iteration <= num_of_nodes:
        e = hq.heappop(q_table)
        if e[2] == iteration:
            selected_nodes.append(e[1])
            curr_spread = calculate_spread(m, selected_nodes)
            iteration += 1
        else:
            margin = -1*(calculate_spread(m, selected_nodes+[e[1]]) - curr_spread)
            iter = iteration
            hq.heappush(q_table, (margin, e[1], iter))
    return selected_nodes


def run():
    m = LTModel(nodes_file, edges_file, status_file, NUM_PLAYERS, player_id-1)

    ns = method(m, nodes_num_per_iter)

    print("Choose: ",ns)
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