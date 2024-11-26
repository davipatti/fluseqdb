import argparse
from pathlib import Path
from fluseqdb import FluSeqDatabase


def main():
    parser = argparse.ArgumentParser(
        "fs_dumpseq_clade", description="Write fasta file comprising a clade."
    )
    parser.add_argument(
        "--segments",
        help="Segment(s) to identify clades in. If multiple segments are passed then the union of "
        "clade members is output. If this argument is omitted, then all segments are used.",
        required=False,
        default=None,
        nargs="+",
    )
    parser.add_argument(
        "--outgroup", help="Tree is rooted on the outgroup.", required=True
    )
    parser.add_argument(
        "--ingroup",
        help="Descendants of the MRCA of the ingroup are included in the output.",
        nargs="+",
    )
    parser.add_argument(
        "--out_dir", help="Directory to write fasta files for each segment."
    )
    parser.add_argument(
        "--db_dir", help="Root directory of a fluseq db. Default='.'.", default="."
    )
    args = parser.parse_args()

    fsdb = FluSeqDatabase(args.db_dir)
    clade = fsdb.clade_members(
        outgroup=args.outgroup, ingroup=args.ingroup, segments=args.segments
    )

    for segment in fsdb.segments:
        fsdb.segments[segment].write_fasta(
            subset=clade, path=Path(args.out_dir, f"{segment}.fasta")
        )
