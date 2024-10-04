import argparse
from fluseqdb import FluSeqDatabase


def main():
    parser = argparse.ArgumentParser(
        "fs_check", description="Check a fluseq db is internally consistent."
    )
    parser.add_argument("--db_dir", help="Root directory of a fluseq db.", default=".")
    args = parser.parse_args()

    FluSeqDatabase(args.db_dir).check()
