# AUTOGENERATED! DO NOT EDIT! File to edit: module6.ipynb (unless otherwise specified).

__all__ = ['Node', 'LinkedList', 'bs']

# Cell
import numpy as np

from .utils import *

# Cell
class Node:
    def __init__(self, value):
        self.value = value
        self.next = None

# Cell
class LinkedList:
    def __init__(self):
        self.head = None
        self.L = 0

    def add(self, value):
        new_node = Node(value)
        head_before = self.head
        self.head = new_node
        new_node.next = head_before
        self.L += 1

    def add_last(self, value):
        new_node = Node(value)
        node = self.head
        while node.next is not None:
            node = node.next
        node.next = new_node
        self.L += 1

    def remove_first(self):
        head = None
        if self.head:
            head = self.head.value
            self.head = self.head.next
            self.L -= 1
        return head

    def __str__(self):
        s = "LinkedList: "
        node = self.head
        while node:
            s += str(node.value) + " ==> "
            node = node.next
        return s

    def __len__(self):
        return self.L

# Cell
def bs(X, val, idx_left, idx_right):
    if idx_left == idx_right:
        if X[idx_left] == val: return idx_left
        else:                  return -1
    else:
        idx_mid = (idx_left + idx_right)//2
        if X[idx_mid] < val:   return bs(X, val, idx_mid+1, idx_right)
        else:                  return bs(X, val, idx_left,  idx_mid)