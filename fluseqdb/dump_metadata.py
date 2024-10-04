import sys
import argparse
from fluseqdb import FluSeqDatabase


def main():
    parser = argparse.ArgumentParser(
        "fsdb_dumpmeta", description="Dump metadata as CSV to stdout"
    )
    parser.add_argument("db_dir", help="Root directory of a fluseq db.")
    args = parser.parse_args()

    FluSeqDatabase(args.db_dir).metadata.sort_values(["source", "isolate_id"]).to_csv(
        sys.stdout
    )
