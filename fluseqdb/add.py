import argparse
import re
from pathlib import Path
from multiprocessing import Pool
from functools import partial

from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord

from .fluseqdb import FluSeqDatabase

REQUIRED_META_DATA_FIELDS = {
    "isolate_id",
    "segment",
}
KNOWN_META_DATA_FIELDS = {
    "isolate_id",
    "dna_insdc",
    "segment",
    "segment_number",
    "isolate_name",
    "collection_date",
    "identifier",
    "dna_accession",
}
SEGMENT_SYNONYMS = {
    "HA_H5": "HA",
    "NA_N1": "NA",
    "NS1": "NS",
}
IGNORE_SEGMENTS = {
    "NS2",  # Subset of NS1
    "PA-X",  # Subset of PA
    "PB1-F2",  # Subset of PB1
    "M2",  # Subset of MP
}


def parse_header(header: str, pattern: str) -> dict[str, str]:
    metadata = re.match(pattern, header).groupdict()

    # Filter out unnecessary metadata
    metadata = {k: v for k, v in metadata.items() if k in KNOWN_META_DATA_FIELDS}

    # Check required metadata is present
    if missing := REQUIRED_META_DATA_FIELDS - set(metadata):
        raise ValueError(f"Missing required metadata: {missing}")

    # Update segment names for synonyms
    if metadata["segment"] in SEGMENT_SYNONYMS:
        current_segment_name = metadata["segment"]
        metadata["segment"] = SEGMENT_SYNONYMS[current_segment_name]

    return metadata


def process_record(record: SeqRecord, fsdb: FluSeqDatabase, header_pattern: str) -> int:
    """
    Process records, adding them to the database if necessary.
    """
    metadata = parse_header(header=record.description, pattern=header_pattern)
    isolate_id = metadata["isolate_id"]
    segment = metadata["segment"]
    sequence_in_db = fsdb.exists(isolate_id=isolate_id, segment=segment)
    if segment not in IGNORE_SEGMENTS and not sequence_in_db:
        fsdb.add(sequence=record.seq, metadata=metadata)


def add(db_dir: str, new_data_dir: str) -> None:
    if not Path(db_dir).is_dir():
        raise ValueError(f"{db_dir} not a directory")

    if not Path(new_data_dir).is_dir():
        raise ValueError(f"{new_data_dir} not a directory")

    fsdb = FluSeqDatabase(db_dir)

    with open(Path(new_data_dir, "header_pattern.txt"), "r") as fobj:
        pattern = "".join(line.strip() for line in fobj.readlines())

    pool = Pool()

    # find fasta file in new_data_dir
    for fasta in Path(new_data_dir).glob("*.fasta"):
        with open(fasta, "r") as fobj:
            records = list(parse(fobj, "fasta"))

        fun = partial(process_record, fsdb=fsdb, header_pattern=pattern)

        pool.map(fun, records)


def main():
    parser = argparse.ArgumentParser(
        "fs_addseqs", description="Add sequences to a fluseq db."
    )
    parser.add_argument(
        "new_data_dir",
        help="Directory containing new data to add. Must contain one or more fasta files and "
        "`header_pattern.txt`",
    )
    parser.add_argument(
        "--db_dir", help="Root directory of a fluseq db. Default='.'.", default="."
    )
    args = parser.parse_args()

    add(db_dir=args.db_dir, new_data_dir=args.new_data_dir)
