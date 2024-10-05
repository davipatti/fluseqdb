from operator import attrgetter
import sys
from string import ascii_lowercase
import random
import os
import subprocess
import json
import logging
from functools import cached_property
from pathlib import Path
from typing import Generator, Optional

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from dendropy import Tree

LOGLEVEL = os.environ.get("LOGLEVEL", "WARNING").upper()
logging.basicConfig(level=LOGLEVEL)


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

    def clade_members(
        self,
        outgroup: str,
        ingroup: list[str],
        include_outgroup: bool = True,
        segments: Optional[list[str]] = None,
    ) -> set[str]:
        """
        Return a set of isolate_ids that are descendants of the MRCA of the ingroup,
        optionally including the outgroup.

        Args:
            outgroup: The isolate_id to root the tree on.
            ingroup: The list of isolate_ids whose MRCA we want to identify descendants of.
            include_outgroup: If True, then the outgroup is included in the set of descendants.
            segments: The segment(s) to identify clades in. If None, then all segments are used.

        Returns:
            A set of isolate_ids that are descendants of the MRCA of the ingroup.
        """
        members = set()
        segments = self.segments if segments is None else segments
        for segment in segments:
            members.update(
                self.segments[segment].clade_members(
                    outgroup=outgroup,
                    ingroup=ingroup,
                    include_outgroup=include_outgroup,
                )
            )
        return members


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

    def write_fasta(
        self, subset: Optional[list[str]] = None, path: Optional[str] = None
    ) -> None:
        """
        Write a fasta file.

        Args:
            subset: Only include isolate_ids in subset in the output.
            path: Path to write the fasta file to disk. If None, then print to stdout.
        """
        subset_intersection = subset & set(record.isolate_id for record in self.records)

        seqs = (
            [
                self.record(isolate_id).sequence
                for isolate_id in sorted(subset_intersection)
            ]
            if subset is not None
            else [
                record.sequence
                for record in sorted(self.records, key=attrgetter("isolate_id"))
            ]
        )

        if path is None:
            write(seqs, sys.stdout, "fasta")
        else:
            with open(path, "w") as fobj:
                write(seqs, fobj, "fasta")

    def tree(self, outgroup: str) -> "dendropy.Tree":
        if not self.record(outgroup):
            raise ValueError(f"{outgroup} not a record")

        uid = "".join(random.choices(ascii_lowercase, k=8))
        tmp_dir = f".{self.segment}_tree_{uid}"
        os.mkdir(tmp_dir)

        ali_path = Path(tmp_dir, "ali.fasta")
        tree_path = Path(tmp_dir, "tree.tre")
        log_path = Path(tmp_dir, "fasttree.log")

        with open(ali_path, "w") as fobj:
            write([record.sequence for record in self.records], fobj, "fasta")

        with (
            open(ali_path, "r") as stdin,
            open(tree_path, "w") as stdout,
            open(log_path, "w") as stderr,
        ):
            subprocess.run(
                ["fasttree", "-nt", "-gtr", "-nosupport"],
                stdin=stdin,
                stdout=stdout,
                stderr=stderr,
            )

        tree = Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
        og = tree.find_node_with_taxon_label(outgroup)
        tree.reroot_at_node(og)
        return tree

    def clade_members(
        self, outgroup: str, ingroup: list[str], include_outgroup: bool = True
    ):
        logging.info(f"building {self.segment} tree...")
        tree = self.tree(outgroup=outgroup)
        clade = extract_mrca_descendants(tree, taxa=ingroup)
        members = set(leaf.taxon.label for leaf in clade.leaf_node_iter())
        if include_outgroup:
            members.add(outgroup)
        return members


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


def extract_mrca_descendants(tree: Tree, taxa: list[str]) -> Tree:
    """
    Extract a subtree comprising all descendants of the MRCA of a set of taxa.
    """
    taxa_in_tree = set(leaf.taxon.label for leaf in tree.leaf_nodes()) & set(taxa)
    mrca = tree.mrca(taxon_labels=list(taxa_in_tree))
    return Tree(seed_node=mrca.extract_subtree())
