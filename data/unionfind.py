
# A fast version of union find using path compression and
# rank-based merging
class UFFast:
    def __init__(self, N):
        self._parent = list(range(N))
        self._rank = [1]*N
    
    def root(self, i):
        """
        Follow parent pointers until reaching a root
        Parameters
        ----------
        i: int
            The starting node 
        
        Returns
        -------
        The root node of i
        """
        p = i
        while self._parent[p] != p:
            p = self._parent[p]
        j = i
        # Path compression
        while self._parent[j] != j:
            j = self._parent[j]
            self._parent[j] = p
        return p

    def get_set_label(self, i):
        """
        Return a number that is the same for every element in
        the set that i is in, and which is unique to that set
        Parameters
        ----------
        i: int
            Element we're looking for
        
        Returns
        -------
        Index of the bubble containing i
        """
        return self.root(i)
    
    def find(self, i, j):
        """
        Return true if i and j are in the same component, or
        false otherwise
        Parameters
        ----------
        i: int
            Index of first element
        j: int
            Index of second element
        """
        return self.root(i) == self.root(j)
    
    def union(self, i, j):
        """
        Merge the two sets containing i and j, or do nothing if they're
        in the same set
        Parameters
        ----------
        i: int
            Index of first element
        j: int
            Index of second element
        """
        root_i = self.root(i)
        root_j = self.root(j)
        if root_i != root_j:
            if self._rank[root_i] < self._rank[root_j]:
                self._parent[root_i] = root_j
                self._rank[root_j] += self._rank[root_i]
            else:
                self._parent[root_j] = root_i
                self._rank[root_i] += self._rank[root_j]