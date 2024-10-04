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

    def __repr__(self):
        return f"FluSeqDatabase(root_dir={self.root_dir})"

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

    def check(self) -> None:
        """
        Check:
            - isolate_ids are consistent in names of sequence files and headers in sequence files.
        """
        for segment in self.segments:
            self.segments[segment].check()


class Record:
    def __init__(self, segment: str, isolate_id: str, fluseqdb: FluSeqDatabase) -> None:
        self.db = fluseqdb
        self.segment = segment
        self.isolate_id = isolate_id

    def __repr__(self):
        return f"Record(segment={self.segment}, isolate_id={self.isolate_id}, fluseqdb={self.db})"

    @cached_property
    def sequence(self) -> SeqRecord:
        path = self.db.path("sequences", self.segment, f"{self.isolate_id}.fasta")
        with open(path) as fobj:
            seq = next(parse(fobj, "fasta"))
            seq._fasta_path = path
            return seq

    @cached_property
    def metadata(self) -> dict:
        path = self.db.path("sequences", self.segment, f"{self.isolate_id}.json")
        with open(path) as fobj:
            return json.load(fobj)

    def check(self) -> None:
        assert (
            self.sequence._fasta_path.stem == self.sequence.id
        ), f"isolate_id in fasta file doesn't match record it contains: {self.sequence._fasta_path}"

        assert (
            self.sequence.id == self.metadata["isolate_id"]
        ), f"isolate_id doesn't match metadata: {self.sequence._fasta_path}"


class SegmentData:

    def __init__(self, fluseqdb: FluSeqDatabase, segment: str) -> None:
        """
        Sequences from a particular segment in a database.
        """
        self.db = fluseqdb
        self.segment = segment

    def __repr__(self):
        return f"SegmentData(fluseqdb={self.db}, segment={self.segment})"

    @property
    def fasta_paths(self):
        return self.db.path("sequences", self.segment).glob("*.fasta")

    @property
    def json_paths(self):
        return self.db.path("sequences", self.segment).glob("*.json")

    def record(self, isolate_id: str) -> Record:
        return Record(self.segment, isolate_id, self.db)

    @property
    def records(self) -> Generator[Record, None, None]:
        for fasta_path in self.fasta_paths:
            isolate_id = fasta_path.stem
            yield self.record(isolate_id)

    @cached_property
    def reference(self) -> SeqRecord:
        """
        A segment's reference sequence.
        """
        with open(self.db.path("reference", f"{self.segment}.fasta")) as fobj:
            return next(parse(fobj, "fasta"))

    @property
    def sequences(self) -> Generator[SeqRecord, None, None]:
        for path in self.fasta_paths:
            with open(path) as fobj:
                seq = next(parse(fobj, "fasta"))
                seq._fasta_path = path
                yield seq

    @cached_property
    def metadata(self) -> pd.DataFrame:
        def generate_rows():
            for path in self.json_paths:
                with open(path) as fobj:
                    yield json.load(fobj)

        return pd.DataFrame(list(generate_rows()))

    def check(self) -> None:
        # Check that records are self-consistent
        for record in self.records:
            record.check()

        # Check no sequence files exist without metadata and vice versa
        fasta_stems = set(path.stem for path in self.fasta_paths)
        json_stems = set(path.stem for path in self.json_paths)
        if missing := json_stems - fasta_stems:
            logging.warn(f"{self.segment} metadata missing sequences: {missing}")
        if missing := fasta_stems - json_stems:
            logging.warn(f"{self.segment} sequences missing metadata: {missing}")


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
    aligner.target_internal_open_gap_score = -5
    aligner.query_internal_open_gap_score = -5

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
