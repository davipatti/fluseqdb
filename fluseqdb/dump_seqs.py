import argparse

from fluseqdb import FluSeqDatabase


def main():
    parser = argparse.ArgumentParser(
        "fs_dumpseq", description="Dump sequences from a fluseq db."
    )
    parser.add_argument(
        "--id_file", help="Optional file of sequence IDs to print, one per line."
    )
    parser.add_argument(
        "--segment", help="Print sequences of this segment.", required=True
    )
    parser.add_argument("--db_dir", help="Root directory of a fluseq db.", default=".")
    args = parser.parse_args()

    if args.id_file is not None:
        with open(args.id_file) as fobj:
            subset = set(line.strip() for line in fobj.readlines())
    else:
        subset = None

    FluSeqDatabase(args.db_dir).segments[args.segment].write_fasta(subset)
