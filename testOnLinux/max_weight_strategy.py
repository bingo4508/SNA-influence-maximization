
import networkx as nx

# N nodes with highest sum of out-egdes influence value
def run(graph, num_of_nodes):
    # write your code here
    nodes_list = graph.nodes(data=True)
    candidate_list = list()

    for n in nodes_list:
        # sum over all influence value of out-edges
        sum_of_out_influence = 0.0
        if n[1]['status'] == 'inactivated':
            for out_edge in graph.out_edges(n[0], data=True):
                sum_of_out_influence += out_edge[2]['influence']
        candidate_list.append((n[0], sum_of_out_influence))

    candidate_list.sort(key = lambda x: x[1], reverse = True)

    # store the selected nodes into return_nodes_list 
    return_num = min(num_of_nodes, len(candidate_list))
    return_nodes_list = list()
    for i in range(0, return_num):
        return_nodes_list.append(candidate_list[i][0])

    return return_nodes_list

