# This Makefile is the template for you to edit. You have to modify it so that we can run your strategy.
# Please note that all parameters could change.
PLAYER_ID=2
NODES_FILE=./networks/egofb_lt_nodes.txt
EDGES_FILE=./networks/egofb_lt_edges.txt
STATUS_FILE=game_status.txt
NODES_NUM_PER_ITER=10
SELECTED_NODES_FILE=selected_nodes.txt
TIME_LIMIT_IN_SEC=60.0

strategy1:
	python3 sample_strategy.py $(PLAYER_ID) $(NODES_FILE) $(EDGES_FILE) $(STATUS_FILE) $(NODES_NUM_PER_ITER) $(SELECTED_NODES_FILE) $(TIME_LIMIT_IN_SEC)

strategy2:
	javac -cp ./jung/* Strategy.java
	java -cp .:./jung/* Strategy $(PLAYER_ID) $(NODES_FILE) $(EDGES_FILE) $(STATUS_FILE) $(NODES_NUM_PER_ITER) $(SELECTED_NODES_FILE) $(TIME_LIMIT_IN_SEC)

strategy3:
	gcc strategy.c -I$(HOME)/.local/include/igraph -L$(HOME)/.local/lib -ligraph -o strategy 
	./strategy $(PLAYER_ID) $(NODES_FILE) $(EDGES_FILE) $(STATUS_FILE) $(NODES_NUM_PER_ITER) $(SELECTED_NODES_FILE) $(TIME_LIMIT_IN_SEC)

clean:
	rm game_status.txt
	rm selected_nodes.txt
