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
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
random.seed(9001)
from random import randint
import statistics

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
    """-	build_kmer_dict /2 qui prend un fichier fastq, une taille k- mer
    et retourne un dictionnaire ayant pour clé le k-mer
    et pour valeur le nombre d’occurrence de ce k-mer"""
    dictionary = {}
    for sequence in read_fastq(fastq_file) :
        for kmer in cut_kmer(sequence, kmer_size):
            if kmer in dictionary:
                dictionary[kmer]+=1
            else :
                dictionary[kmer]=1
    return dictionary


def build_graph(kmer_dict):
    """La fonction build_graph /4 prendra en entrée un dictionnaire de k-mer
    et créera l’arbre de k-mers préfixes et suffixes décrit précédemment.
    Les arcs auront pour paramètre obligatoire un poids nommé “weight”.

Validez à l’aide de pytest le test de construction, testez votre implémentation de pylint et n’oubliez pas de commiter et pusher.
"""
    G=nx.DiGraph()
    for kmer in kmer_dict:
        kmer1=""
        kmer2=""
        for i in range(len(kmer)-1):
            kmer1 = kmer1+kmer[i]
            kmer2 = kmer2+kmer[i+1]
        G.add_edge(kmer1, kmer2, weight=kmer_dict[kmer])
    return G


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass


def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass


def path_average_weight(graph, path):
    pass


def solve_bubble(graph, ancestor_node, descendant_node):
    pass

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
    for start in starting_nodes :
        for end in ending_nodes :
            for contig in nx.all_simple_paths(graph, start, end, cutoff=None) :
                string_contig = contig[0][:-1] + "".join([sommet[0] for sommet in contig[1:-1]]) + contig[-1]
                contigs.append((string_contig, len(string_contig)))
    return contigs
    # (début, node milieu.., node fin, nombre de nodes)

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """-	save_contigs /2 qui prend une liste de tuple (contig, taille du contig)
    et un nom de fichier de sortie et écrit un fichier de sortie contenant les contigs selon le format fasta
    (retour chariot tous les 80 caractères) à l’aide de la fonction fill:"""



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
    print(myfile_path)

    graph = build_graph(build_kmer_dict(myfile_path, 3))
    draw_graph(graph, "graphimg_file")
    print(get_starting_nodes(graph))
    get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))

    print("++++"*30)
    graph = nx.DiGraph()
    graph.add_edges_from(
        [("TC", "CA"), ("AC", "CA"), ("CA", "AG"), ("AG", "GC"), ("GC", "CG"), ("CG", "GA"), ("GA", "AT"),
         ("GA", "AA")])
    contig_list = get_contigs(graph, ["TC", "AC"], ["AT", "AA"])
    print("contig_list print :", contig_list)
    # main()

