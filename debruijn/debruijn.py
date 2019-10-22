# -*- coding: utf-8 -*-

"""
@author : Zygnematophyce
Master II BIB - 2019 2020
Projet cours: Réalisation d'un assembleur.
"""

import argparse
import networkx as nx


def arguments():
    """ Méthode qui définie les arguments avec argparse."""
    parser = argparse.ArgumentParser(description='Debruij.py')
    parser.add_argument("-i",
                        help="fichier fastq single end",
                        type=str)
    parser.add_argument("-k",
                        help="-k taille des kmer (optionnel - default 21)",
                        type=int)
    parser.add_argument("-o",
                        help="-o fichier config",
                        type=str)
    args = parser.parse_args()

    return args.i, args.k, args.o

def cut_kmer(sequence, length_kmer):
    """ Méthode qui retourne la séquence  de k-mer."""
    for i in range(0, len(sequence) - length_kmer + 1, 1):
        yield sequence[i:i+length_kmer]

def read_fastq(fastq_path):
    """ Méthode qui retourne une séquence du fastq."""
    with open(fastq_path, "r") as fastq_file:
        for line in fastq_file.readlines():
            if line[0] == "A" or line[0] == "T" or line[0] == "G" or line[0] == "C":
                yield line.strip("\n")

def build_kmer_dict(fastq_path, length_kmer):
    """
    Méthode qui retourne un dictionnaire avec comme clès le k-mer et
    valeur le nombre d'occurence de ce k-mer.
    """
    k_mer_dict = dict()

    for seq in read_fastq(fastq_path):
        for kmer in cut_kmer(seq, length_kmer):
            if kmer in k_mer_dict.keys():
                k_mer_dict[kmer] += 1
            else:
                k_mer_dict[kmer] = 1
    return k_mer_dict

def build_graph(dictionnary_kmer):
    """
    Méthode qui prend en entrée un dicitonnaire de k-mer
    et créé un arbre de k-mer.
    """
    print("build_graph")

    # Créer un graph dirigé.
    graph = nx.DiGraph()
    for key in dictionnary_kmer:
        graph.add_edge(u_of_edge=key[:-1],
                       v_of_edge=key[1:],
                       weight=dictionnary_kmer[key])
    return graph

def get_starting_nodes(graph_network):
    """
    Méthode qui prend en entrée un graphe et retourne une liste
    de noeuds d'entrée.
    """
    return list(graph_network.nodes)

def get_sink_nodes(graph_network):
    """
    Méthode qui prend en entrée un graphe et retourne une liste
    de noeuds de sortie.
    """
    print("get sink nodes")

def get_contigs(input_graph_network,
                output_graph_network):
    """
    prend un graphe, une liste de noeuds d’entrée et une liste de sortie et
    retourne une liste de tuple(contig, taille du contig)
    """
    print("get contigs")

def save_contigs(graph_tuple, output_name):
    """
    qui prend un tuple (contig, taille du contig) et un nom de fichier de
    sortie et écrit un fichier de sortie contenant les contigs selon le format.
    """
    print("save contigs")

def std():
    pass

def path_average_weight():
    pass

def remove_paths():
    pass

def select_best_path():
    pass

def solve_bubble():
    pass

def simplify_bubbles():
    pass

def solve_entry_tips():
    pass

def solve_out_tips():
    pass

# main
if __name__ == "__main__":
    print("Debruij")

    # Definition des variables avec les arguments.
    FASTQ_FILE, KMER_LENGTH, OUTPUT_FILE = arguments()

    """
    Retourne un dictionnaire du motif du k_mere avec son nombre d'occurence
    contenu dans des fastq.
    """
    KMER_DICT = build_kmer_dict(FASTQ_FILE, KMER_LENGTH)

    # Construction du graph avec le dictionnaire.
    network = build_graph(KMER_DICT)

    # Construction d'une liste de noeuds a partir d'un graph.
    nodes = get_starting_nodes(network)

    # Afficher tous les noeuds du grap.
    print(nodes)
