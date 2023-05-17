# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2023/4/8 15:47
# @Author     : zzhen
# @File       : node.py
# @Software   : PyCharm
# @Description: 
# @Copyright  : Copyright (c) 2023 by zzhen, All Rights Reserved.

from functools import lru_cache
import logging

logger = logging.getLogger('easy353.DBG')


class Node:
    # use __slots__ to reduce memory usage
    __slots__ = ('value', 'occur', 'parents', 'children', "len")

    def __init__(self, value: int, occurrence: int, k_size: int, parents: set = None, children: dict = None):
        # value is binary representation of the kmer
        self.value = value
        self.occur = occurrence
        self.len = k_size
        # type:set
        self.parents = parents if parents else set()
        # type:dict
        self.children = children if children else {}

    def __repr__(self):
        return f'Node({self.value},{self.occur},{self.len},{len(self.parents)},{len(self.children)})'

    def __str__(self):
        return f'Node({self.value},{self.occur},{self.len},{len(self.parents)},{len(self.children)})'

    # if the value and the length of two nodes are the same, they are equal
    def __eq__(self, other):
        return self.value == other.value and self.len == other.len

    # rewrite hash function
    # if the value and the length of two nodes are the same, they have the same hash value
    def __hash__(self):
        return hash((self.value, self.len))

    # add node to the set of children
    # add self to the set of parents of the out node
    def add_child(self, node):
        self.children[node.value] = node
        node.parents.add(self)

    # delete node from the set of children
    def del_child(self, node):
        self.children.pop(node.value)
        node.parents.remove(self)

    # Use LRU Cache to speed up the method
    @lru_cache(maxsize=1000)
    def get_children(self, value):
        return self.children.get(value)


# stupid function: when backtracking, the path becomes shorter and shorter
def path_effectiveness(path_a: list, path_b: list, k_size: int, threshold: float = 0.1):
    def path_stats(path: list, size: int):
        length = sum(elem.len - size + 2 for elem in path) + size - 2
        occur = sum(elem.occur for elem in path)
        return length, occur

    path_a_len, path_a_occur = path_stats(path_a, k_size)
    path_b_len, path_b_occur = path_stats(path_b, k_size)
    # if one path has longer length but lower occurrence,return the other path
    if (path_a_len >= path_b_len) != (path_a_occur >= path_b_occur):
        return path_a if path_a_occur >= path_b_occur else path_b
    # calculate the len_effectiveness of the length
    len_effectiveness = 2 * abs(path_a_len - path_b_len) / (path_a_len + path_b_len)
    # if the len_effectiveness of the two paths is greater than the threshold, return the longer path
    if len_effectiveness >= threshold:
        return path_a if path_a_len > path_b_len else path_b
    # if the len_effectiveness of the two paths is less than the threshold, return the path with higher avg occurrence
    return path_a if path_a_occur / path_a_len > path_b_occur / path_b_len else path_b


class DeBruijnGraph:
    # use __slots__ to reduce memory usage
    __slots__ = ('nodes', 'k_size')

    def __init__(self, k_size):
        self.nodes = {}
        self.k_size = k_size

    def __repr__(self):
        graph_str = "De Bruijn Graph:\n"
        for node in self.nodes.values():
            graph_str += node.__repr__()
        return graph_str

    # build the graph from a list of kmers
    # split each kmer into prefix and suffix
    # node1: prefix, node2: suffix and add edge between them
    def build_de_bruijn_graph(self, read_kmer_dict):
        for bin_kmer, occur in read_kmer_dict.items():
            # the prefix of a kmer is the first k-1 bases
            prefix = bin_kmer >> 2
            # the suffix of a kmer is the last k-1 bases
            suffix = bin_kmer & ((1 << (2 * (self.k_size - 1))) - 1)
            # if prefix or suffix not in the graph, add them
            # otherwise, update their occurrence
            if prefix not in self.nodes:
                self.nodes[prefix] = Node(prefix, occur, self.k_size - 1)
            else:
                self.nodes[prefix].occur += occur
            if suffix not in self.nodes:
                self.nodes[suffix] = Node(suffix, occur, self.k_size - 1)
            else:
                self.nodes[suffix].occur += occur
            # add edge between prefix and suffix
            self.nodes[prefix].add_child(self.nodes[suffix])

    # remove low occurrence nodes
    def simple_de_bruijn_graph(self):
        for node in self.nodes.values():
            # remove low occurrence children nodes
            if len(node.children) <= 1:
                continue
            # the threshold is 20% of the occurrence of the most frequent child node
            occur_threshold = max(node.children.values(), key=lambda x: x.occur).occur * 0.2
            children_to_remove = [child for child in node.children.values() if child.occur < occur_threshold]
            for child in children_to_remove:
                node.del_child(child)

    # compress the dBG: if the node has only one parent and child node, merge them
    def dbg_compress_nodes(self):
        nodes_to_add, nodes_to_remove = {}, set()
        # Iterate all nodes in the dBG
        for node in self.nodes.values():
            # if the node has been visited, skip it
            if node in nodes_to_remove:
                continue
            # if the node has only one parent and child node, merge them
            if len(node.children) == 1 and len(node.parents) == 1:
                # get the parent and child node
                p_node = next(iter(node.parents))
                c_node = next(iter(node.children.values()))
                if len(p_node.children) == 1 and len(c_node.parents) == 1:
                    # the value of the new node is the concatenation of the three nodes
                    tmp_node_value = (p_node.value << (
                            2 * ((node.len + c_node.len) - 2 * (self.k_size - 2)))) | (
                                             node.value << (2 * (c_node.len - self.k_size + 2))) | c_node.value
                    # the occurrence of the new node is the max occurrence of the three nodes
                    tmp_node_occur = sum([p_node.occur, node.occur, c_node.occur])
                    tmp_node_len = p_node.len + node.len + c_node.len - 2 * (self.k_size - 2)
                    # If the new value is not in the compressed_nodes, add it
                    # In fact, this condition is always true
                    if tmp_node_value not in nodes_to_add:
                        # Add the nodes to be the nodes_to_remove set
                        # compressed_node may have the same value with c_node
                        nodes_to_remove.update({p_node, node, c_node})
                        # Create the new node
                        # The parents and children of the new node are the parents and children of the merged nodes
                        compressed_node = Node(tmp_node_value, tmp_node_occur, tmp_node_len,
                                               parents=p_node.parents, children=c_node.children)
                        nodes_to_add[tmp_node_value] = compressed_node
                        # Update the children and parents of the merged nodes
                        for parent in p_node.parents:
                            parent.children.pop(p_node.value)
                            parent.add_child(compressed_node)
                        for child in c_node.children.values():
                            child.parents.remove(c_node)
                            child.parents.add(compressed_node)
        # Update the compressed nodes in the dBG
        self.nodes.update(nodes_to_add)
        # Remove the nodes in the nodes_to_remove set
        for node in nodes_to_remove:
            # add the if condition to avoid the error: KeyError
            # because the node with the same value may have been removed in the previous iteration
            # eg.Node(803690428455,3,20) Node(803690428455,6,22)
            if node.value in self.nodes:
                del self.nodes[node.value]

    @staticmethod
    # the function using DFS and backtracking to traverse the graph
    def dfs_backtrack_traversal(init_node: Node, max_iterations: int = 1000, forward: bool = True):
        visited_nodes, stack, best_path, cur_path, node_distance = [init_node], [], [], [], 0
        while True:
            # the direction of the path: forward or backward
            # if forward, the next nodes are the children nodes
            # if backward, the next nodes are the parents nodes
            next_nodes = sorted(visited_nodes[-1].children.values() if forward else visited_nodes[-1].parents,
                                key=lambda elem: elem.occur, reverse=True)
            # filter out visited nodes
            next_nodes = [node for node in next_nodes if node not in visited_nodes]
            # if there are more than one next_nodes, add them to the stack for backtracking
            # record the node_distance and reset it to 0
            if len(next_nodes) >= 2:
                stack.append((next_nodes[1:], node_distance))
                node_distance = 0
            # if there is no next nodes, turn back to the previous node and update the best path
            if not next_nodes:
                max_iterations -= 1
                # if the len of the current path is longer than the best path, update the best path
                if sum([elem.len for elem in cur_path]) > sum([elem.len for elem in best_path]):
                    best_path = cur_path.copy()
                # pop the nodes which have been visited in this iteration
                for _ in range(node_distance):
                    visited_nodes.pop(), cur_path.pop()
                if stack:
                    # 如果直接pop出上一个多父节点或者多子节点的节点，每一个节点的子节点或父节点最多只被访问两个，这样会导致路径不完整
                    next_nodes, node_distance = stack.pop()
                # if there is no nodes in the stack_list, break the loop
                else:
                    break
            visited_nodes.append(next_nodes[0])
            cur_path.append(next_nodes[0])
            node_distance += 1
            # when the max_iterations is reached, break the loop
            if not max_iterations:
                break
        return best_path if forward else best_path[::-1]
