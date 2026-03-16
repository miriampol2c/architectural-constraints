import argparse
from CLEAN.utils import *


def eval_parse():
    # only argument passed is the fasta file name to infer
    # located in ./data/[args.fasta_data].fasta
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--fasta_data', type=str)
    args = parser.parse_args()
    return args


def main():
    args = eval_parse()
    # esm embedding are taken care of
    ensure_dirs("data/esm_data")
    ensure_dirs("data/pretrained")
    # fasta file must be in this case alocated in data/
    retrive_esm1b_embedding(args.fasta_data)

if __name__ == '__main__':
    main()

