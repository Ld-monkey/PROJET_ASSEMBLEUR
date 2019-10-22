# -*- coding: utf-8 -*-

"""
@author : Zygnematophyce
Master II BIB - 2019 2020
Projet cours: Réalisation d'un assembleur.
"""

import argparse


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
            if line[0] == "A" or line[0] == "T" or line[0] == "G" or line[0] =="C":
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
            else :
                k_mer_dict[kmer] = 0
    return k_mer_dict

# main
if __name__ == "__main__":
    print("Debruij")

    # Definition des variables avec les arguments.
    FASTQ_FILE, KMER_LENGTH, OUTPUT_FILE = arguments()

    """
    Retourne un dictionnaire du motif du k_mere avec son nombre d'occurence
    contenu dans des fastq.
    """
    kmer_dict = build_kmer_dict(FASTQ_FILE, KMER_LENGTH)
    print(kmer_dict)
