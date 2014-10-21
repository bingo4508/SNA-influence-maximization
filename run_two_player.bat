set nodes_file=./networks/egofb_lt_nodes.txt
set edges_file=./networks/egofb_lt_edges.txt
set p1_strategy_id=1
set p2_strategy_id=1
set iter=10
set nodes_num_per_iter=10
set first_time_limit=300
set time_limit=60

C:\Python34\python.exe main_two_player_test.py %nodes_file% %edges_file% %p1_strategy_id% %p2_strategy_id% %iter% %nodes_num_per_iter% %first_time_limit% %time_limit%
