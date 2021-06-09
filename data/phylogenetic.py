import numpy as np
import matplotlib.pyplot as plt
import json
from unionfind import *

class TreeNode(object):
    """
    Attributes
    ----------
    left: TreeNode
        Left child
    right: TreeNode
        Right child
    sim: int
        Phylogenetic similarity between two children if
        internal node, or max similarity between any two
        nodes if leaf node
    key: string
        Name of species if leaf node
    """

    def __init__(self, param):
        """
        Initialize a tree node

        Parameters
        ----------
        param: string or int
            If string, this is a leaf node with the name of a species
            If int, then this is an internal node with a particular
            height in the phylogenetic tree
        """
        self.sim = 0
        self.key = None
        if type(param) is str:
            self.key = param
        else:
            self.sim = param
        self.left = None
        self.right = None
    
    def __str__(self):
        """
        String is key if key is not None, or blank otherwise
        """
        ret = ""
        if self.key:
            ret = "{}".format(self.key)
        return ret

    def compute_y_coords(self, maxsim=[0], y=[0]):
        """
        Recursively compute y coordinate of nodes via an inorder
        traversal, while computing the maximum phylogenetic 
        similarity as a side effect

        Parameters
        ----------
        maxsim: list of [int]
            Maximum similarity
        y: list of [int]
            Current y coordinate
        """
        if self.left:
            self.left.compute_y_coords(maxsim, y)
        maxsim[0] = max(maxsim[0], self.sim)
        self.y = y[0]
        y[0] += 1
        if self.right:
            self.right.compute_y_coords(maxsim, y)

    def compute_x_coords(self, maxsim):
        """
        Recursively compute and store the x coordinates
        of all nodes.  If the nodes are internal, then the
        x coordinate is the phylogenetic similarity.
        If the node is a leaf node, then the x coordinate
        is the maximum phylogenetic similarity among all
        internal nodes

        Parameters
        ----------
        maxsim: int
            Maximum phylogenetic similarity across all nodes
        """
        if self.left:
            self.left.compute_x_coords(maxsim)
        if self.right:
            self.right.compute_x_coords(maxsim)
        if self.key:
            self.x = maxsim
        else:
            self.x = self.sim
            

    def draw(self):
        """
        Recursively draw phylogenetic tree.  Assumes that the
        x and y coordinates have been precomputed
        """
        x1, y1 = self.x, self.y
        # Draw a dot
        plt.scatter(x1, y1, 50, 'k')
        # Draw some text indicating what the key is
        plt.text(x1+10, y1, "{}".format(self))
        if self.left:
            # Draw a line segment from my node to this left child
            x2, y2 = self.left.x, self.left.y
            plt.plot([x1, x2], [y1, y2])
            self.left.draw()
        if self.right:
            # Draw a line segment from my node to this right child
            x2, y2 = self.right.x, self.right.y
            plt.plot([x1, x2], [y1, y2])
            self.right.draw()

class PhyloTree(object):
    def __init__(self):
        self.root = None

    def draw(self, threshold=None):
        """
        Draw the phylogenetic tree from the bottom up

        Parameters
        ----------
        threshold: int
            If specified, draw a vertical line showing a similarity
            threshold for clustering
        """
        if self.root:
            maxsim = [0]
            self.root.compute_y_coords(maxsim)
            self.root.compute_x_coords(maxsim[0])
            self.root.draw()
            ax = plt.gca()
            xlim = ax.get_xlim()
            ax.set_xlim([xlim[0], xlim[1]+200])
            if threshold:
                ylim = ax.get_ylim()
                plt.plot([threshold, threshold], [ylim[0], ylim[1]], 'k', linestyle='--', linewidth=3)
                plt.title("Similarity Threshold = {}".format(threshold))
            ax.set_yticks([])
            plt.xlabel("Needleman-Wunsch Similarity")
            plt.tight_layout()


def load_blosum(filename):
    """
    Load in a BLOSUM scoring matrix for Needleman-Wunsch

    Parameters
    ----------
    filename: string
        Path to BLOSUM file
    
    Returns
    -------
    A dictionary of {string: int}
        Key is string, value is score for that particular 
        matching/substitution/deletion
    """
    fin = open(filename)
    lines = [l for l in fin.readlines() if l[0] != "#"]
    fin.close()
    symbols = lines[0].split()
    X = [[int(x) for x in l.split()] for l in lines[1::]]
    X = np.array(X, dtype=int)
    N = X.shape[0]
    costs = {}
    for i in range(N-1):
        for j in range(i, N):
            c = X[i, j]
            if j == N-1:
                costs[symbols[i]] = c
            else:
                costs[symbols[i]+symbols[j]] = c
                costs[symbols[j]+symbols[i]] = c
    return costs
