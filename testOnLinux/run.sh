#!/bin/bash 
nodes_file=./networks/egofb_lt_nodes.txt
edges_file=./networks/egofb_lt_edges.txt 
strategy_id=1
iter=10
nodes_num_per_iter=10
first_time_limit=300
time_limit=60
python3 main.py $nodes_file $edges_file $strategy_id $iter $nodes_num_per_iter $first_time_limit $time_limit
