#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import random
import os
import argparse
import sys
import statistics
import networkx as nx
import numpy as np
import tqdm

random.seed(9001)

__author__ = "Aurore Wilbrink"
__copyright__ = "CY Tech - Cergy"
__credits__ = ["Aurore Wilbrink"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aurore Wilbrink"
__email__ = "wilbrinkau@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
    "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract from the file Sequences
      :Parameters:
          fastq_file: file .fq
    Returns: A generator of sequences
    """
    with open(fastq_file, 'rt') as myfile:
        sentences = myfile.readlines()
        for sentence in sentences:
            sentence = sentence[:-1]
            if (
                    sentence[0] == "T"
                    or sentence[0] == "G"
                    or sentence[0] == "A"
                    or sentence[0] == "C"
            ):
                yield sentence


def cut_kmer(read, kmer_size):
    """ Generator of k-mer
      :Parameters:
          read : sequence extracted with read_fastq
          kmer_size : size of our k-mer
    Returns: K-mer
    """
    for i in range(len(read) - kmer_size + 1):
        kmer = read[i:i + kmer_size]
        yield kmer


def build_kmer_dict(fastq_file, kmer_size):
    """ Creation of a kmer dictionary
      :Parameters:
          fastq_file : a file fastq
          kmer_size : size of kmers
    Returns: a dictionary of k-mer {"key" = number of occurancy}
    """
    dictionary = {}
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, kmer_size):
            if kmer in dictionary:
                dictionary[kmer] += 1
            else:
                dictionary[kmer] = 1
    return dictionary


def build_graph(kmer_dict):
    """ Creation of graph
      :Parameters:
          kmer_dict : dictionary of kmer
    Returns: a graph according to the kmer's dictionary
    """

    graph = nx.DiGraph()
    for kmer in kmer_dict:
        kmer1 = ""
        kmer2 = ""
        for i in range(len(kmer) - 1):
            kmer1 = kmer1 + kmer[i]
            kmer2 = kmer2 + kmer[i + 1]
        graph.add_edge(kmer1, kmer2, weight=kmer_dict[kmer])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ Function in order to remove paths
      :Parameters:
          graph : DiGraph
          path_list : list of path
          delete_entry_node : boolean in order to delete or not the entry_node
          delete_sink_node : boolean in order to delete (or not) the sinking_node
    Returns: a new graph without the unwilling paths
    """
    for path in path_list:
        for i in range(len(path)):
            if not (delete_entry_node) and i == 0:
                continue
            if not (delete_sink_node) and i == len(path) - 1:
                continue
            try:
                graph.remove_node(path[i])
            except nx.exception.NetworkXError:
                pass
    return graph


def std(data):
    """ Calculate standard deviation
      :Parameters:
          data : a number
    Returns: Standard deviation of this number
    """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """
    Select the best paths and clean the graph from its unecessary paths
      :Parameters:
        graph : Digraph
        path_list : paths that we want to compare (list)
        path_length : length of each path (list (int))
        weight_avg_list : weight average of each path (list (int))
        delete_entry_node : delete entry node of the path (bool)
        delete_sink_node : delete sink node of the path (bool)
      Returns:
        graph : Digraph cleaned
    """

    index = []
    max_weight = max(weight_avg_list)

    for i in range(len(path_list)):
        if weight_avg_list[i] < max_weight:
            graph = remove_paths(graph, [path_list[i]], delete_entry_node, delete_sink_node)
            index.append(i)
    index = sorted(index, reverse=True)
    for indice in index:
        path_list.pop(indice)
        path_length.pop(indice)
        weight_avg_list.pop(indice)

    if len(path_list) == 1:
        return graph

    max_length = max(path_length)
    index = []
    for i in range(len(path_list)):
        if path_length[i] != max_length:
            graph = remove_paths(graph, [path_list[i]], delete_entry_node, delete_sink_node)
            index.append(i)
    index = sorted(index, reverse=True)
    for indice in index:
        path_list.pop(indice)
        path_length.pop(indice)
        weight_avg_list.pop(indice)

    if len(path_list) == 1:
        return graph

    choice = random.randint(0, len(path_list))

    for i in range(len(path_list)):
        if i != choice:
            graph = remove_paths(graph, [path_list[i]], delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    """
    Calculate the average weight of a specific path
    :param graph: Digraph
    :param path: path of nodes
    :return: average weight
    """
    weights = []
    for i in range(len(path) - 1):
        sommet_depart = path[i]
        sommet_arrive = path[i + 1]
        weights.append(graph[sommet_depart][sommet_arrive]["weight"])
    return np.mean(weights)


def solve_bubble(graph, ancestor_node, descendant_node):
    """
    Remove bubble from the graph
    :param graph: Digraph
    :param ancestor_node: node
    :param descendant_node: node
    :return: graph without bubble
    """
    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    if paths==[]:
        return graph
    graph = select_best_path(
        graph,
        paths,
        [len(path) for path in paths],
        [path_average_weight(graph, path) for path in paths]
    )
    return graph


def simplify_bubbles(graph):
    """
    Removes every bubbles from a graph
    :param graph: Digraph
    :return: Digraph
    """

    for start in tqdm.tqdm(get_starting_nodes(graph)):
        for end in tqdm.tqdm(get_sink_nodes(graph)):
            graph = solve_bubble(graph, start, end)
    return graph


def solve_entry_tips(graph, starting_nodes):
    """
    Removes unecessary starting nodes in a graph
    :param graph: Digraph
    :param ending_nodes: list of starting nodes
    :return: graph cleaned of unecessary starting nodes
    """
    return solve_out_tips(graph.reverse(), starting_nodes).reverse()


def solve_out_tips(graph, ending_nodes):
    """
    Removes unecessary ending nodes in a graph
    :param graph: Digraph
    :param ending_nodes: list of ending nodes
    :return: graph cleaned of unecessary ending nodes
    """
    for i in tqdm.tqdm(range(len(ending_nodes))):
        end_1 = ending_nodes[i]
        for j, end_2 in enumerate(ending_nodes):
            if i >= j:
                continue
            if not (end_1 in graph) or not (end_2 in graph):
                continue
            ancestor_node = nx.algorithms.lowest_common_ancestors.lowest_common_ancestor(
                graph,
                end_1,
                end_2)
            if ancestor_node == None:
                continue
            paths = list(nx.all_simple_paths(graph, ancestor_node, end_1)) + list(
                nx.all_simple_paths(graph, ancestor_node, end_2))
            if paths == []:
                continue
            graph = select_best_path(
                graph,
                paths,
                [len(path) for path in paths],
                [path_average_weight(graph, path) for path in paths],
                delete_sink_node=True
            )
    return graph


def get_starting_nodes(graph):
    """
    Get the starting nodes
    :param graph: Digraph
    :return: list of starting_nodes
    """
    starting_nodes = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return np.array(starting_nodes)


def get_sink_nodes(graph):
    """
    get the ending nodes
    :param graph: DiGraph
    :return: list of ending_nodes
    """
    ending_nodes = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)
    return np.array(ending_nodes)


def get_contigs(graph, starting_nodes, ending_nodes):
    """
    Return all the contigs into a list of tuple
    :param graph: DiGraph
    :param starting_nodes: list of starting_nodes
    :param ending_nodes: list of ending nodes
    :return: list of tuple(contig, size_contig)
    """
    contigs = []
    for start in starting_nodes:
        for end in ending_nodes:
            for contig in nx.all_simple_paths(graph, start, end, cutoff=None):
                contig = [str(contig_i) for contig_i in contig]
                string_contig = contig[0][:-1] \
                                + "".join([sommet[0] for sommet in contig[1:-1]]) + contig[-1]
                contigs.append((string_contig, len(string_contig)))
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i + width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """
    Save contigs into a file (fasta format)
    :param contigs_list: list of tuples(contig, size_contig)
    :param output_file: name of the output file
    :return: a file output_file.fq
    """
    file = open(output_file, "w")
    i = 0
    for sequence, lenght in contigs_list:
        line = ">contig_{0} len={1}\n{2}\n".format(str(i), str(lenght), fill(sequence))
        file.write(line)
        i += 1
    file.close()


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = vars(get_arguments())

    print("Lecture du fichier et création du graphe")
    graph = build_graph(build_kmer_dict(args['fastq_file'], args['kmer_size']))
    print("Résolution des bulles")
    graph = simplify_bubbles(graph)
    print("Résolution des pointes d’entrée et de sortie")
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    print("Ecriture du/des contigs")
    save_contigs(
        get_contigs(
            graph,
            get_starting_nodes(graph),
            get_sink_nodes(graph)
        ),
        os.path.join(args['output_file'])
    )


if __name__ == '__main__':
    main()
