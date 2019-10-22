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

    # Affiche les argument
    if args.i:
        print(args.i)
    if args.k:
        print(args.k)
    if args.o:
        print(args.o)

    return args.i, args.k, args.o



# main
if __name__ == "__main__":
    print("Debruij")

    # Definition des variables avec les arguments.
    FASTQ_FILE, KMER_LENGTH, OUTPUT_FILE = arguments()
