import json
import logging
from functools import cached_property
from pathlib import Path
from typing import Generator

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

SEGMENTS = "PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"
SEGMENT_SPECIFIC_METADATA = [
    "dna_insdc",
    "segment",
    "segment_number",
    "dna_accession",
    "identifier",
]


class FluSeqDatabase:

    def __init__(self, root_dir: str) -> None:
        """
        A database of curated influenza sequences.

        Args:
            root_dir: The root directory of the database.
        """
        self.root_dir = root_dir
        self.segments = {segment: SegmentData(self, segment) for segment in SEGMENTS}

    def path(self, *args) -> Path:
        """
        Construct a path relative to the database's root directory.
        """
        return Path(self.root_dir, *args)

    def exists(self, isolate_id: str, segment: str) -> bool:
        """
        Check if isolate_id exists for segment.
        """
        return self.path("sequences", segment, f"{isolate_id}.fasta").exists()

    def add(self, sequence: str | SeqRecord, metadata: dict[str, str]) -> None:
        """
        Add a sequence to the database.
        """
        segment = metadata["segment"]
        isolate_id = metadata["isolate_id"]
        try:
            reference_seq = self.segments[segment].reference
        except KeyError:
            logging.warning(f"unrecognised segment: {metadata}")
            return

        try:
            aligned = align_to_reference(
                reference_seq=reference_seq, input_seq=sequence
            )
        except ValueError:
            logging.warning(
                f"internal gaps in reference when aligning {segment} for {isolate_id}"
            )
            return

        with open(self.path("sequences", segment, f"{isolate_id}.fasta"), "w") as fobj:
            write(
                SeqRecord(seq=Seq(aligned), id=isolate_id, description=""),
                handle=fobj,
                format="fasta",
            )

        with open(self.path("sequences", segment, f"{isolate_id}.json"), "w") as fobj:
            fobj.write(json.dumps(metadata))

    @cached_property
    def metadata(self) -> pd.DataFrame:
        return pd.concat(
            self.segments[segment].metadata.drop(columns=SEGMENT_SPECIFIC_METADATA)
            for segment in self.segments
        ).drop_duplicates()


class SegmentData:

    def __init__(self, fluseqdb: FluSeqDatabase, segment: str) -> None:
        """
        Sequences from a particular segment in a database.
        """
        self.db = fluseqdb
        self.segment = segment

    @cached_property
    def reference(self) -> SeqRecord:
        """
        A segment's reference sequence.
        """
        with open(self.db.path("reference", f"{self.segment}.fasta")) as fobj:
            return next(parse(fobj, "fasta"))

    def load_seqs(self) -> Generator[SeqRecord, None, None]:
        for path in self.db.path("sequences", self.segment).glob("*.fasta"):
            with open(path) as fobj:
                yield next(parse(fobj, "fasta"))

    @cached_property
    def metadata(self) -> pd.DataFrame:
        def generate_rows():
            for path in self.db.path("sequences", self.segment).glob("*.json"):
                with open(path) as fobj:
                    yield json.load(fobj)

        return pd.DataFrame(list(generate_rows()))


def align_to_reference(reference_seq: str, input_seq: str) -> str:
    """
    Align an input sequence to a reference. Returns the aligned input sequence trimmed to the region
    of the reference.

    Args:
        reference_seq (str):
        input_seq (str):

    Raises:
        ValueError: If internal gaps are introduced in the reference during alignment.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"

    # Bigger penalty to open gaps
    aligner.target_internal_open_gap_score = -1e2
    aligner.query_internal_open_gap_score = -1e2

    # Perform the alignment
    alignments = aligner.align(reference_seq, input_seq)
    best_alignment = alignments[0]

    # Extract the aligned sequences
    # Biopython doesn't provide a nice way of extracting the aligned sequences. This is the
    # 'cleanest' I could come up with.
    aligned_ref, _, aligned_input = best_alignment._format_unicode().strip().split("\n")

    # Check for internal gaps in the reference sequence
    first_non_gap, last_non_gap = first_and_last_non_gap(aligned_ref)
    if "-" in aligned_ref[first_non_gap:last_non_gap]:
        msg = "Internal gaps introduced in reference sequence"
        logging.warning(f"{msg}:\n{best_alignment.format()}")
        raise ValueError(msg)

    # Trim the input sequence to match the reference length
    # Get the positions corresponding to the non-gap parts of the reference
    trimmed_input_seq = aligned_input[first_non_gap : last_non_gap + 1]

    return trimmed_input_seq


def first_and_last_non_gap(sequence: str) -> tuple[int, int]:
    """
    Returns the indices of the first and last non-gap ('-') characters in a sequence.
    If all characters are gaps, returns None for both indices.

    Args:
        sequence (str): The input sequence containing gaps ('-').

    Returns:
        tuple: A tuple (first_non_gap_index, last_non_gap_index).
    """
    for i, char in enumerate(sequence):
        if char != "-":
            first_non_gap = i
            break

    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != "-":
            last_non_gap = i
            break

    return first_non_gap, last_non_gap
