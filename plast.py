import argparse
import os



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    #add argument -i:
    parser.add_argument(
        '-i',
        type=str,
        help='input sequence'
    )
    parser.add_argument(
        '-db',
        type=str,
        help='database file (fasta format)',
        default = "tRNAs.fasta"
    )
    parser.add_argument(
        "-E",
        type = int
        )
    parser.add_argument(
        "-ss",
        type = int
        )
    parser.add_argument(
        "-seed",
        type = str
        )
