from functools import cached_property
from operator import attrgetter
from pathlib import Path
from string import ascii_lowercase
from typing import Generator, Optional
import json
import logging
import os
import random
import subprocess
import sys

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord
from dendropy import Tree
import pandas as pd


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
        Initialize a FluSeqDatabase from a root directory.

        The root directory should have a subdirectory "sequences" with subdirectories
        for each segment (PB2, PB1, PA, HA, NP, NA, MP, NS). The root directory should also contain
        a subdirectory "reference" that contains a fasta file for each segment containing a
        reference sequence.

        Args:
            root_dir: The root directory of the database.
        """
        if not (
            Path(root_dir, "reference").exists()
            and Path(root_dir, "sequences").exists()
        ):
            raise ValueError(
                "root directory must contain 'reference' and 'sequences' subdirectories."
            )

        self.root_dir = root_dir
        self.segments = {segment: SegmentData(self, segment) for segment in SEGMENTS}

    def __repr__(self):
        return f"FluSeqDatabase(root_dir={self.root_dir})"

    def path(self, *args) -> Path:
        """
        Return a Path object for the given subdirectory/file in the database root
        directory.

        Args:
            *args: The path components to join to the root directory.

        Returns:
            A Path object for the given path.
        """
        return Path(self.root_dir, *args)

    def exists(self, isolate_id: str, segment: str) -> bool:
        """
        Check if a sequence exists in the database.

        Args:
            isolate_id: The isolate_id of the sequence to check.
            segment: The segment of the sequence to check.

        Returns:
            True if the sequence exists, otherwise False.
        """
        return self.path("sequences", segment, f"{isolate_id}.fasta").exists()

    def add(self, sequence: str | SeqRecord, metadata: dict[str, str]) -> None:
        """
        Add a sequence to the database.

        Args:
            sequence: The sequence to add as either a str or SeqRecord.
            metadata: A dictionary of metadata for the sequence.
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
        """
        A DataFrame of all metadata in the database.

        This DataFrame does not include the columns 'dna_insdc', 'segment', 'segment_number',
        'dna_accession', 'identifier' since these are specific to a particular segment.
        """
        return pd.concat(
            self.segments[segment].metadata[
                ["isolate_id", "isolate_name", "collection_date"]
            ]
            for segment in self.segments
        ).drop_duplicates()

    def check(self) -> None:
        """
        Check all sequences in the database for consistency.

        This method calls the check method of all segment-specific SequenceData objects.
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
        """
        A Record is a particular sequence and its metadata in the database.

        Args:
            segment: The segment of the sequence.
            isolate_id: The isolate_id of the sequence.
            fluseqdb: The database the sequence is from.
        """
        self.db = fluseqdb
        self.segment = segment
        self.isolate_id = isolate_id

    def __repr__(self):
        return f"Record(segment={self.segment}, isolate_id={self.isolate_id}, fluseqdb={self.db})"

    @cached_property
    def sequence(self) -> SeqRecord:
        """
        The sequence as a SeqRecord.

        The _fasta_path attribute is set to the path of the fasta file this record
        was read from.
        """
        path = self.db.path("sequences", self.segment, f"{self.isolate_id}.fasta")
        with open(path) as fobj:
            seq = next(parse(fobj, "fasta"))
            seq._fasta_path = path
            return seq

    @cached_property
    def metadata(self) -> dict:
        """
        The metadata associated with this record as a dictionary.

        The metadata is loaded from disk the first time it is accessed and then
        cached for subsequent accesses.
        """
        path = self.db.path("sequences", self.segment, f"{self.isolate_id}.json")
        with open(path) as fobj:
            return json.load(fobj)

    def check(self) -> None:
        """
        Check the record is self-consistent.

        Checks that the isolate_id in the fasta filename matches the isolate_id in
        the fasta record and that the isolate_id in the fasta record matches the
        isolate_id in the metadata.

        Raises an AssertionError if either check fails.
        """
        assert (
            self.sequence._fasta_path.stem == self.sequence.id
        ), f"isolate_id in fasta file doesn't match record it contains: {self.sequence._fasta_path}"

        assert (
            self.sequence.id == self.metadata["isolate_id"]
        ), f"isolate_id doesn't match metadata: {self.sequence._fasta_path}"


class SegmentData:

    def __init__(self, fluseqdb: FluSeqDatabase, segment: str) -> None:
        """
        A collection of sequences and their metadata for a particular segment.

        Args:
            fluseqdb: The database the sequences are from.
            segment: The segment of the sequences.
        """
        self.db = fluseqdb
        self.segment = segment

    def __repr__(self):
        return f"SegmentData(fluseqdb={self.db}, segment={self.segment})"

    @property
    def fasta_paths(self):
        """
        A list of paths to the fasta files for the sequences in this segment.

        Returns:
            Generator of Path objects to the fasta files.
        """
        return self.db.path("sequences", self.segment).glob("*.fasta")

    @property
    def json_paths(self):
        """
        A list of paths to the json files for the sequences in this segment.

        Returns:
            Generator of Path objects to the json files.
        """
        return self.db.path("sequences", self.segment).glob("*.json")

    def record(self, isolate_id: str) -> Record:
        """
        Get a Record for a particular isolate_id in this segment.

        Args:
            isolate_id: The isolate_id of the sequence to get.

        Returns:
            A Record object for the sequence.
        """
        return Record(self.segment, isolate_id, self.db)

    @property
    def records(self) -> Generator[Record, None, None]:
        """
        A generator of Record objects for all sequences in this segment.

        Yields:
            Generator of Record objects for the sequences in this segment.
        """
        for fasta_path in self.fasta_paths:
            isolate_id = fasta_path.stem
            yield self.record(isolate_id)

    @cached_property
    def reference(self) -> SeqRecord:
        """
        The reference sequence for this segment.

        Returns:
            A SeqRecord object for the reference sequence for this segment.
        """
        with open(self.db.path("reference", f"{self.segment}.fasta")) as fobj:
            return next(parse(fobj, "fasta"))

    @property
    def sequences(self) -> Generator[SeqRecord, None, None]:
        """
        A generator of SeqRecord objects for all sequences in this segment.

        The sequences are read from disk on demand and the _fasta_path attribute is set
        to the path of the fasta file the sequence was read from.

        Yields:
            Generator of SeqRecord objects for the sequences in this segment.
        """
        for path in self.fasta_paths:
            with open(path) as fobj:
                seq = next(parse(fobj, "fasta"))
                seq._fasta_path = path
                yield seq

    @cached_property
    def metadata(self) -> pd.DataFrame:
        """
        A DataFrame of all metadata for all sequences in this segment.

        The DataFrame is generated from the .json files in the segment directory
        and is cached for subsequent accesses.

        Returns:
            A DataFrame of all metadata for all sequences in this segment.
        """

        def generate_rows():
            for path in self.json_paths:
                with open(path) as fobj:
                    yield json.load(fobj)

        return pd.DataFrame(list(generate_rows()))

    def check(self) -> None:
        # Check that records are self-consistent
        """
        Check all sequences in this segment for consistency.

        Checks that all sequences are self-consistent (i.e. the isolate_id in the fasta
        filename matches the isolate_id in the fasta record and that the isolate_id in
        the fasta record matches the isolate_id in the metadata) and that no sequence
        files exist without metadata and vice versa.
        """
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
        Write all sequences in this segment to a fasta file.

        If a subset of isolate_ids is specified, only those sequences are written.
        If a path is specified, the sequences are written to that file. Otherwise,
        the sequences are written to stdout.

        Args:
            subset: A list of isolate_ids to write. Default is None, which causes
                all sequences to be written.
            path: The path to write the sequences to. Default is None, which causes
                the sequences to be written to stdout.
        """
        if subset is None:
            seqs = [
                record.sequence
                for record in sorted(self.records, key=attrgetter("isolate_id"))
            ]
        else:
            intersection = subset & set(record.isolate_id for record in self.records)
            seqs = [
                self.record(isolate_id).sequence for isolate_id in sorted(intersection)
            ]

        if path is None:
            write(seqs, sys.stdout, "fasta")
        else:
            with open(path, "w") as fobj:
                write(seqs, fobj, "fasta")

    def tree(self, outgroup: str) -> Tree:
        """
        Compute a tree for all sequences in this segment, rooted at the given outgroup.

        Args:
            outgroup: The isolate_id of the sequence to root the tree on.

        Returns:
            A dendropy.Tree object representing the tree.
        """
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
        tree.reroot_at_edge(og.edge)
        return tree

    def clade_members(
        self, outgroup: str, ingroup: list[str], include_outgroup: bool = True
    ):
        """
        Return a set of isolate_ids that are descendants of the MRCA of the ingroup,
        optionally including the outgroup.

        Args:
            outgroup: The isolate_id to root the tree on.
            ingroup: The list of isolate_ids whose MRCA we want to identify descendants of.
            include_outgroup: If True, then the outgroup is included in the set of descendants.

        Returns:
            A set of isolate_ids that are descendants of the MRCA of the ingroup.
        """
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
    Extracts the subtree of a tree that descends from the MRCA of a given set of taxa.

    Args:
        tree (Tree): The tree to extract the subtree from.
        taxa (list[str]): The set of taxa to find the MRCA of.

    Returns:
        Tree: The subtree of `tree` that descends from the MRCA of `taxa`.

    Notes:
        The returned subtree will include the MRCA itself as the root node.
    """
    taxa_in_tree = set(leaf.taxon.label for leaf in tree.leaf_nodes()) & set(taxa)
    mrca = tree.mrca(taxon_labels=list(taxa_in_tree))
    return Tree(seed_node=mrca.extract_subtree())
