__author__ = 'bingo4508'


def write_nodes_list(nodes_list, selected_nodes_file):
    with open(selected_nodes_file, 'w') as f:
        for n in nodes_list:
            f.write('%d ' % n)