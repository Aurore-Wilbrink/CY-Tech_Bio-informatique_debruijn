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
from pathlib import Path
import statistics
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
random.seed(9001)



__author__ = "Aurore Wilbrink"
__copyright__ = "CY Tech - Cergy"
__credits__ = ["Aurore Wilbrink"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aurore Wilbrink"
__email__ = "wilbrinkau@eisti.eu"
__status__ = "Developpement"

debug = True

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
            if sentence[0] == "T" or sentence[0] == "G" or sentence[0] == "A" or sentence[0] == "C" :
                yield sentence


def cut_kmer(read, kmer_size):
    """ Generator of k-mer
      :Parameters:
          read : sequence extracted with read_fastq
          kmer_size : size of our k-mer
    Returns: K-mer
    """
    for i in range(len(read)-kmer_size + 1):
        kmer = read[i:i+kmer_size]
        yield kmer


def build_kmer_dict(fastq_file, kmer_size):
    """ Creation of a kmer dictionary
      :Parameters:
          fastq_file : a file fastq
          kmer_size : size of kmers
    Returns: a dictionary of k-mer {"key" = number of occurancy}
    """
    dictionary = {}
    for sequence in read_fastq(fastq_file) :
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
        for i in range(len(kmer)-1):
            kmer1 = kmer1+kmer[i]
            kmer2 = kmer2+kmer[i+1]
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
            if not(delete_entry_node) and i == 0:
                continue
            if not(delete_sink_node) and i == len(path) - 1:
                continue
            graph.remove_node(path[i])
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
    """qui prend un graphe, une liste de chemin, une liste donnant la longueur de chaque chemin,
    une liste donnant le poids moyen de chaque chemin, delete_entry_node pour indiquer
    si les noeuds d’entrée seront supprimés
    et delete_sink_node pour indiquer si les noeuds de sortie seront supprimés
    et retourne un graphe nettoyé des chemins indésirables.
    Par défaut, delete_entry_node et delete_sink_node seront ici à False."""

    index = []
    max_weight = max(weight_avg_list)

    print(len(path_list))
    for i in range(len(path_list)):
        print("i = ", i)
        print("weight_avg_list[i] : ", weight_avg_list[i])
        print("max_weight : ", max_weight)
        if weight_avg_list[i] != max_weight :
            print("[path_list[i]] : ", [path_list[i]])
            graph = remove_paths(graph, [path_list[i]], delete_entry_node, delete_sink_node)
            index.append(i)
    for indice in index :
        path_list.pop(indice)
        print("path_list après : ", path_list)
        print("path_length avant : ", path_length)
        path_length.pop(indice)
        print("path_length après : ", path_length)
        print("weight_avg_list avant : ", weight_avg_list)
        weight_avg_list.pop(indice)
        print("weight_avg_list après : ", weight_avg_list)

    if len(path_list) == 1:
        return graph

    max_length = max(path_length)
    index = []
    print("Avant length : len path list : ", len(path_list))
    for i in range(len(path_list)):
        if path_length[i] != max_length:
            print("[path_list[i]] : ", [path_list[i]])
            graph = remove_paths(graph, [path_list[i]], delete_entry_node, delete_sink_node)
            index.append(i)
    for indice in index:
        print("path_list avant : ", path_list)
        path_list.pop(indice)
        print("path_list après : ", path_list)
        print("path_length avant : ", path_length)
        path_length.pop(indice)
        print("path_length après : ", path_length)
        print("weight_avg_list avant : ", weight_avg_list)
        weight_avg_list.pop(indice)
        print("weight_avg_list après : ", weight_avg_list)

    if len(path_list) == 1:
        return graph

    choice = random.randint(0,len(path_list))

    for i in range(len(path_list)):
        if i != choice :
            print("[path_list[i]] : ", [path_list[i]])
            graph = remove_paths(graph, [path_list[i]], delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    """qui prend un graphe et un chemin et qui retourne un poids moyen."""
    weights = []
    for i in range(len(path)-1):
        sommet_depart = path[i]
        sommet_arrive = path[i+1]
        weights.append(graph[sommet_depart][sommet_arrive]["weight"])
    return np.mean(weights)


def solve_bubble(graph, ancestor_node, descendant_node):
    """qui prend un graphe, un noeud ancêtre, un noeud descendant
    et retourne un graph nettoyé de la bulle se trouvant entre ces deux noeuds
    en utilisant les fonctions précédemment développée."""
    for noeud in graph : print(noeud)
    graph = remove_paths(graph, get_contigs(graph, [ancestor_node], [descendant_node]), delete_entry_node=False, delete_sink_node=False)
    print("new")
    for noeud in graph: print(noeud)
    return graph

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    """prend en entrée un graphe et retourne une liste de noeuds d’entrée """
    # générateur tous les successeurs ou prefecesseur, on les compte
    starting_nodes = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return np.array(starting_nodes)


def get_sink_nodes(graph):
    """prend en entrée un graphe et retourne une liste de noeuds de sortie """
    ending_nodes = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)
    return np.array(ending_nodes)


def get_contigs(graph, starting_nodes, ending_nodes):
    """prend un graphe, une liste de noeuds d’entrée et une liste de sortie
    et retourne une liste de tuple(contig, taille du contig) """
    contigs = []
    if debug:
        print("starting node :", starting_nodes)
        print("ending_nodes : ", ending_nodes)
    for start in starting_nodes:
        for end in ending_nodes:
            for contig in nx.all_simple_paths(graph, start, end, cutoff=None):
                if debug: print("contig = ", contig)
                string_contig = contig[0][:-1] + "".join([sommet[0] for sommet in contig[1:-1]]) + contig[-1]
                contigs.append((string_contig, len(string_contig)))
    if debug: print("contigs = ", contigs)
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """-	save_contigs /2 qui prend une liste de tuple (contig, taille du contig)
    et un nom de fichier de sortie et écrit un fichier de sortie contenant les contigs selon le format fasta
    (retour chariot tous les 80 caractères) à l’aide de la fonction fill:"""
    # print(contigs_list)

    f = open(output_file, "w")
    i = 0
    for sequence, lenght in contigs_list:
        line = ">contig_{0} len={1}\n{2}\n".format(str(i), str(lenght), fill(sequence))
        f.write(line)
        i += 1
    f.close()

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    path = Path(__file__).resolve().parent
    myfile_path = Path("..") / "eva71_two_reads.fq"
    graph = build_graph(build_kmer_dict(myfile_path, 4))

    save_contigs(get_contigs(graph,get_starting_nodes(graph), get_sink_nodes(graph)), "output_file.txt")
    # main()

    graph = nx.DiGraph()
    graph.add_edges_from(
        [("TCA", "CAA"), ("ACA", "CAA"), ("CAA", "AGA"), ("AGA", "GCA"), ("GCA", "CGA"), ("CGA", "GAA"), ("GAA", "ATA"),
         ("GAA", "AAA")])
    contig_list = get_contigs(graph, ["TCA", "ACA"], ["ATA", "AAA"])


    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from([(1, 2, 10), (3, 2, 10), (2, 4, 15),
                                     (4, 5, 15), (2, 10, 10), (10, 5, 10),
                                     (2, 8, 3), (8, 9, 3), (9, 5, 3),
                                     (5, 6, 10), (5, 7, 10)])
    graph_1 = solve_bubble(graph_1, 2, 5)
