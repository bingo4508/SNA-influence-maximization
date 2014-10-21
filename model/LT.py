__author__ = 'bingo4508'

import networkx as nx
import copy

class LTModel():
    def __init__(self, nodes_file, edges_file, status_file, player_num, my_id):
        self.my_id = my_id
        self.player_num = player_num
        self.G = None
        self.newly_active_nodes_list = []
        self.active_nodes_list = []
        self.selected_nodes = []
        self.rival_selected_nodes = [[]]*self.player_num
        self.snapshot = {}
        self.load_graph(nodes_file, edges_file, status_file)

    # load graph from files
    def load_graph(self, nodes_file, edges_file, status_file):
        self.G = nx.DiGraph()
        # Nodes
        with open(nodes_file, 'r') as f:
            f.readline()
            for l in f:
                l = l.rstrip().split(' ')
                self.G.add_node(int(l[0]), status='inactivated', threshold=0, owner=None,
                                influenced=[0.0 for i in range(0, self.player_num)], test=False)
                self.G.node[int(l[0])]['threshold'] = float(l[1])

        # Edges
        with open(edges_file, 'r') as f:
            f.readline()
            f.readline()
            for l in f:
                l = l.rstrip().split(' ')
                self.G.add_edge(int(l[0]), int(l[1]), influence=float(l[2]))

        # Update Status
        with open(status_file, 'r') as f:
            owner = 0
            queue_num = self.player_num+self.my_id
            queue = [[]]*queue_num
            quene_index = 0
            for line in f:
                nodes = [int(n) for n in line.rstrip().split(' ')]
                quene_index = (quene_index+1) % queue_num
                queue[quene_index] = nodes
                self.rival_selected_nodes[owner] = nodes
                for n in nodes:
                    self.G.node[int(n)]['owner'] = owner
                    self.G.node[int(n)]['status'] = 'activated'
                    self.active_nodes_list.append(n)
                    owner = (owner+1) % self.player_num
            # Set current newly activated nodes
            for e in queue:
                for n in e:
                    self.G.node[int(n)]['status'] = 'newly_activated'
                    self.newly_active_nodes_list.append(n)
                    self.active_nodes_list.remove(n)

    def select_nodes(self, nodes_list, test=False):
        self.selected_nodes = nodes_list
        for n in nodes_list:
            if self.G.node[n]['status'] == 'inactivated':
                # Save original state in snapshot
                if test:
                    self.snapshot[n] = copy.deepcopy(self.G.node[n])

                self.G.node[n]['status'] = 'newly_activated'
                self.G.node[n]['owner'] = self.my_id
                self.newly_active_nodes_list.append(n)

    def propagate(self, test=False):
        while True:
            copy_newly_active_nodes = self.newly_active_nodes_list[:]
            self.newly_active_nodes_list = []
            if copy_newly_active_nodes:
                for n in copy_newly_active_nodes:
                    # Save original state in snapshot
                    if test:
                        if n not in self.snapshot:
                            self.snapshot[n] = copy.deepcopy(self.G.node[n])

                    self.G.node[n]['status'] = 'activated'
                    self.active_nodes_list.append(n)
                    owner = self.G.node[n]['owner']
                    for e in self.G.successors(n):
                        if self.G.node[e]['status'] == 'inactivated':
                            # Save original state in snapshot
                            if test:
                                if e not in self.snapshot:
                                    self.snapshot[e] = copy.deepcopy(self.G.node[e])

                            self.G.node[e]['influenced'][owner] += self.G.edge[n][e]['influence']

                            if self.G.node[e]['influenced'][owner] >= self.G.node[e]['threshold']:
                                self.G.node[e]['status'] = 'newly_activated'
                                self.newly_active_nodes_list.append(e)
                                # Determine the owner
                                success = True
                                for p in range(0, self.player_num):
                                    if self.G.node[e]['influenced'][owner] < self.G.node[e]['influenced'][p]:
                                        success = False
                                if success:
                                    self.G.node[e]['owner'] = owner
            else:
                break

    def reset(self):
        anl = self.active_nodes_list

        for k, v in self.snapshot.items():
            if self.G.node[k]['status'] == 'activated':
                anl.remove(k)
            self.G.node[k] = v

        self.snapshot = {}

    def check(self, g):
        for n in self.G.nodes():
            if self.G.node[n]['status'] != g.node[n]['status'] or \
                self.G.node[n]['owner'] != g.node[n]['owner'] or \
                self.G.node[n]['influenced'][0] != g.node[n]['influenced'][0] or \
                self.G.node[n]['influenced'][1] != g.node[n]['influenced'][1]:
                print("ERROR!!!")
                print("Current: ", n, self.G.node[n])
                print("Original: ", n, g.node[n])
                break

    def get_activated_nodes(self):
        c = [0]*self.player_num
        for n in self.active_nodes_list:
            c[self.G.node[n]['owner']] += 1
        return c

    # return a copy of graph
    def get_copy_graph(self):
        return self.G.copy()