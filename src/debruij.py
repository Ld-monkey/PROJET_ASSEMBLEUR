import argparse


# main
if __name__ == "__main__":
    print("Debruij")

    """
    -i fichier fastq single end
    -k taille des kmer (optionnel - default 21)
    -r Reference genome (optionnel)
    -o fichier config
    """

    # Création d'un parser
    parser = argparse.ArgumentParser(description = 'Debruij.py')
    parser.add_argument("-i", "--fichier_fastq", help = "fichier fastq single end",
                        type="")

    args = parser.parse_args()

    if args.fichier_fastq:
        print(args.fichier_fastq)

