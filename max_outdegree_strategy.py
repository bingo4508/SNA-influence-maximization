
import networkx as nx

# N nodes with highest out-degree
def run(graph, num_of_nodes):
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
