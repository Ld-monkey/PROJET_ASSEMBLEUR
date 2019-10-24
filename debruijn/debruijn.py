# -*- coding: utf-8 -*-

"""
@author : Zygnematophyce
Master II BIB - 2019 2020
Projet cours: Réalisation d'un assembleur.
"""

import argparse
import networkx as nx
import os
import statistics
import random

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

def build_kmer_dict(fastq_path, length_kmer=21):
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
    input_node_list = list()
    for input_nodes in graph_network.nodes:
        if len(list(graph_network.predecessors(input_nodes))) == 0:
            input_node_list.append(input_nodes)
    return input_node_list

def get_sink_nodes(graph_network):
    """
    Méthode qui prend en entrée un graphe et retourne une liste
    de noeuds de sortie.
    """
    output_node_list = list()
    for output_nodes in graph_network.nodes:
        if len(list(graph_network.successors(output_nodes))) == 0:
            output_node_list.append(output_nodes)
    return output_node_list

def get_contigs(network_graph,
                input_graph_network,
                output_graph_network):
    """
    Méthode qui retourne une liste de tuple(contig, taille du contig)
    """
    contigs = []
    for noeud_depart in input_graph_network:
        for noeud_fin in output_graph_network:
            for path in nx.all_simple_paths(network_graph, source=noeud_depart, target=noeud_fin):
                prep_contig = path
                contig_ecrit = []
                contig_ecrit.append(prep_contig[0])
                for i in range(1, len(prep_contig)):
                    contig_ecrit.append(prep_contig[i][-1:])
                contig_ecrit = "".join(contig_ecrit)
                contigs.append((contig_ecrit, len(contig_ecrit)))
    return contigs


def fill(text, width=80):
    """split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(graph_tuple, output_name):
    """
    qui prend un tuple (contig, taille du contig) et un nom de fichier de
    sortie et écrit un fichier de sortie contenant les contigs selon le format.
    """
    numero = 0
    print("save contigs")
    with open(output_name, "w") as output_file:
        for contigs in graph_tuple:
            output_file.write(">contig_{} len={}\n".format(numero, contigs[1]))
            output_file.write("{}\n".format(fill(contigs[0])))
            numero += 1

def std(list_value):
    """ Méthode qui retourne l'écart type. """
    return statistics.stdev(list_value)

def path_average_weight(network_graph, path_list):
    """
    Méthode qui prend un graphe et un chemin et qui retourne un poids
    moyen.
    """
    list_poids = list()
    for u,v,e in network_graph.subgraph(path_list).edges(data=True):
        list_poids.append(e['weight'])
    return statistics.mean(list_poids)

def remove_paths(network_graph, paths, delete_entry_node, delete_sink_node):
    """
    Qui prend un graphe et une liste de chemin, delete_entry_node pour
    indiquer si les noeuds d’entrée seront supprimés et delete_sink_node
    pour indiquer si es noeuds de sortie seront supprimés et retourne un
    graphe nettoyé des chemins indésirables.
    """
    graph = network_graph
    for i in range(len(paths)):
        graph.remove_nodes_from(paths[i][1:-1])
        if delete_entry_node == True:
            graph.remove_node(paths[i][0])
        if delete_sink_node == True:
            graph.remove_node(paths[i][-1])
    return graph


def select_best_path(network_graph, paths, delete_entry_node, delete_sink_node):
    """
    Qui prend un graphe et une liste de chemin, delete_entry_node pour
    indiquer si les noeuds d’entrée seront supprimés et delete_sink_node
    pour indiquer si les noeuds de sortie seront supprimés et retourne un
    graphe nettoyé des chemins indésirables.
    """
    """
    par 3 critères:
    - Un chemin est plus fréquent
    - Un chemin est plus long
    - Le hasard, vous imposerez une seed à 9001
    """
    pass

def solve_bubble():
    """
    Qui prend un graphe, un noeud ancêtre, un noeud descendant et
    retourne un graph nettoyé de la bulle se trouvant entre ces
    deux noeuds.
    """

def simplify_bubbles():
    """
    Qui prend un graphe et retourne un graphe sans bulle
    """

def solve_entry_tips():
    """
    Qui prend un graphe et une liste de noeuds d’entrée et retourne
    graphe sans chemin d’entrée indésirable
    """

def solve_out_tips():
    """
    Qui prend un graphe et une liste de noeuds de sortie et retourne graphe
    sans chemin de sortie indésirable
    """

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

    # Construction d'une liste de noeuds d'entrés a partir d'un graph.
    all_input_nodes = get_starting_nodes(network)

    # Construction d'une liste de noeuds de sortie a partir d'un graph.
    all_output_nodes = get_sink_nodes(network)

    # Afficher tous les noeuds du grap.
    print(all_input_nodes)
    print(all_output_nodes)

    # Afficher le graph.
    #nx.draw(network)
    #plt.show()
    all_contigs = list()
    all_contigs = get_contigs(network,
                              all_input_nodes,
                              all_output_nodes)
    print(all_contigs)
    save_contigs(all_contigs, OUTPUT_FILE)

    #path_average_weight(network, )

