import multiprocessing as np
import subprocess as sp
import time
import sys
import networkx as nx

import DiffusionModel

import max_weight_strategy as s1
import max_outdegree_strategy as s2


'''
This is the main procedure of hw1 (the arena of influence maximization game).
You must install python3 and install NetworkX on python3 
If you have any question, please contact TAs of SNA courses.
'''

# read the list of selected nodes from file
def read_nodes_from_file(filename):
    nodes_list = list()
    with open(filename, 'r') as f:
        line = f.readline()
        entry = line.strip().split()
        for e in entry:
            nodes_list.append(int(e))
    return nodes_list

# update the game status to the target file
def update_status_file(filename, nodes_list):
    with open(filename, 'a') as f:
        for n in nodes_list:
            print(n, end=' ', file=f)
        print('',file=f)

# clear the file content
def clear_file(filename):
    with open(filename, 'w') as f:
        pass

# main function 
if __name__ == '__main__':
    if len(sys.argv) < 9:
        print('Usage:', sys.argv[0], 'nodes_file edges_file player1_strategy_id player2_strategy_id max_iter nodes_num_per_iter first_time_limit(secs) time_limit(secs)', file=sys.stderr)
        exit(-1)
    
    # arguments
    nodes_file = sys.argv[1]
    edges_file = sys.argv[2]
    p1_strategy_id = int(sys.argv[3])
    p2_strategy_id = int(sys.argv[4])
    max_iter = int(sys.argv[5])
    nodes_num_per_iter = int(sys.argv[6])
    first_time_limit = float(sys.argv[7])
    time_limit = float(sys.argv[8])
    status_file = 'game_status.txt'
    selected_nodes_file = 'selected_nodes.txt'
    assert max_iter >= 0
    assert time_limit >= 0.0

    # clear the game status file
    clear_file(status_file)
    clear_file(selected_nodes_file)

    # read in graph
    print('Read in graph ...', end='', file = sys.stderr)
    model = DiffusionModel.MultiPlayerLTModel(nodes_file, 
            edges_file, player_num = 2)
    print('Done', file=sys.stderr)

    # main loop
    print('Influence Maxization Game Starts !',file=sys.stderr)
    for i in range(0, max_iter):
        print('\n\n===============Iteration %d===============' % (i), file=sys.stderr)
        # initialize a round 
        model.init_round()
        if i == 0:
            time_limit_this_round = first_time_limit
        else:
            time_limit_this_round = time_limit
		
        ## Player 1 ##
        print('Player 1 is selecting nodes ...',end='', file=sys.stderr)
        try:
            sp.call(['java','-jar','SNA_HW1.jar','1',nodes_file,edges_file,status_file,
			str(nodes_num_per_iter),selected_nodes_file,str(time_limit_this_round)],
			timeout=time_limit_this_round)
        except sp.TimeoutExpired:
            print(' TimeLimitedExceed ... ', end='', file=sys.stderr)
        except:
            print(' RunTimeError ... ', end='', file=sys.stderr)
        
        selected_nodes_list = read_nodes_from_file(selected_nodes_file)
        selected_nodes_list = model.select_nodes(selected_nodes_list, player_id = 0)

        # Output the selected nodes of player 1
        update_status_file(status_file, selected_nodes_list)
        print('Done', file=sys.stderr)
        print('Player 1 successfully selected the following nodes:', selected_nodes_list, file=sys.stderr)


        ## Player 2: your strategy ##
        # copy_graph = model.get_copy_graph()
        # selected_nodes_list = s1.run(copy_graph, nodes_num_per_iter)
        # selected_nodes_list = model.select_nodes(selected_nodes_list, player_id = 1)
        print('Player 2 is selecting nodes ...',end='', file=sys.stderr)
        # try:
            # sp.call(['java','-jar','SNA_HW1.jar','2',nodes_file,edges_file,status_file,
			# str(nodes_num_per_iter),selected_nodes_file,str(time_limit_this_round)],
			# timeout=time_limit_this_round)
        # except sp.TimeoutExpired:
            # print(' TimeLimitedExceed ... ', end='', file=sys.stderr)
        # except:
            # print(' RunTimeError ... ', end='', file=sys.stderr)
        
        selected_nodes_list = read_nodes_from_file(selected_nodes_file)
        selected_nodes_list = model.select_nodes(selected_nodes_list, player_id = 1)

        # Output the selected nodes of player 2
        update_status_file(status_file, selected_nodes_list)
        print('Done', file=sys.stderr)
        print('Player 2 successfully selected the following nodes:', selected_nodes_list, file=sys.stderr)

        # run linear threshold once
        print('Propagate once ...', end='',file=sys.stderr)
        model.propagate()
        update_status_file(status_file, model.get_activated_nodes(player_id=0))
        update_status_file(status_file, model.get_activated_nodes(player_id=1))
        print('Done\n', file=sys.stderr)

        if model.no_nodes_left():
            print('No inactivated node left. The game ends.', file=sys.stderr)
            break
           
        model.print_result()

