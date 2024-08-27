import numpy as np

def RemoveFrom(l, xs):
    for x in xs:
        try:
            l.remove(x)
        except:
            raise Exception(x)

class Network:
    def __init__(self, edges, dists = None, symmetry = False):
        """
        edge dict: 
            key: index of node, int
            val: indices of linked nodes, list of ints
        weight dict:
            key: index of node, int
            val: weights of linked nodes, array of floats
        neighbours dict:
            key: index of node, int
            val: list of lists of ints
        """
        self.nodes = list(edges.keys())
        self.size = len(self.nodes)
        self.edges = edges
        if symmetry:
            for key in self.edges.keys():
                for node in self.edges[key]:
                    if key not in self.edges[node]:
                        self.edges[node].append(key)
        if not (dists is None):
            self.dists = dists
            
        else:
            # init equal weights
            self.dists = {key: np.array([1 for _ in value]) for key, value in edges.items()}
        self.neighbours, self.neighbour_weights = {},{}
        for key, value in edges.items():
            neighbours, neighbour_weights = self.ComputeNeighbours(key)
            self.neighbours[key] = neighbours
            self.neighbour_weights[key] = neighbour_weights
        self.adj_mat = self.AdjMatrix()
        self.w_mats = self.WeightMatrix()
            
    def AdjMatrix(self):
        """
        return the adjacency matrix
        """
        mat = np.zeros((self.size,self.size))
        for i in range(self.size):
            mat[([i for j in self.edges[i]],[j for j in self.edges[i]])] = 1
        return mat
        
        
    def ComputeNeighbours(self, node):
        neighbours = [self.edges[node]]
        neighbour_weights = [1/self.dists[node]/sum(1/self.dists[node])]
        nodes_left = self.nodes[:] # nodes that are not yet connected to this node
        RemoveFrom(nodes_left,neighbours[-1]+[node])
        prev_dists = self.dists[node]
        while True:
            stage_neighbours = []
            stage_dists = {}
            for i in range(len(neighbours[-1])):
                # for each stage_r-1 neighbour, find its neighbours
                prev_neighbour = neighbours[-1][i]
                for j in range(len(self.edges[prev_neighbour])):
                    stage_neighbour = self.edges[prev_neighbour][j]
                    if (stage_neighbour in nodes_left) & (stage_neighbour!=node):
                    # if the stage-r neighbour is not in previous stage(still left unconnected) and is not self
                        if not (stage_neighbour in stage_neighbours):
                        # if stage neighbour not added to stage-r neighbourhood, add to neighbourhood 
                        # and add dist(node to stage neighbour) = min_dist(node to prev neighbour)+dist(prev neighbour to stage neighbour)
                            stage_neighbours.append(stage_neighbour)
                            stage_dists[stage_neighbour] = [prev_dists[i]+self.dists[prev_neighbour][j]]
                        else:
                        # add dist(...) = ....
                            stage_dists[stage_neighbour].append(prev_dists[i]+self.dists[prev_neighbour][j])
            if len(stage_neighbours) == 0:
                # if no new neighbours detected, end and return
                return neighbours, neighbour_weights
            # for each stage-r neighbour, find the min_dist(node to stage neighbour)
            stage_dists = np.array([min(stage_dists[stage_neighbour]) for stage_neighbour in stage_neighbours])
            RemoveFrom(nodes_left, stage_neighbours) # update nodes not connected
            neighbours.append(stage_neighbours)
            neighbour_weights.append(1/stage_dists/sum(1/stage_dists))
            prev_dists = stage_dists

    def WeightMatrix(self):
        max_stage = max([len(self.neighbours[node]) for node in self.nodes])
        w_mats = [np.identity(self.size)]
        for r in range(max_stage):
            mat = np.zeros((self.size,self.size))
            for i in range(self.size):
                if len(self.neighbours[i]) > r:
                    mat[([i for j in self.neighbours[i][r]],[j for j in self.neighbours[i][r]])] = 1
            w_mats.append(np.divide(mat, np.sum(mat,axis=0), out=np.zeros_like(mat), where=np.sum(mat,axis=0)!=0))
        return np.array(w_mats)

n8_14=Network({0:[1,3,4],1:[0,6],2:[3,4],3:[2,4],4:[0,2],5:[0,1,4],6:[3],7:[4,6,2]},symmetry=True)
n8_7=Network({0:[1],1:[0],2:[4],3:[4],4:[2],5:[0,1],6:[3],7:[6]},symmetry=True)
n8_7_2=Network({0:[6,7],1:[3,4],2:[5],3:[1,5],4:[1,6],5:[2,3],6:[0,4],7:[0]},symmetry=True)