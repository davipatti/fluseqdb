# fluseqdb

A simple influenza sequence database.

## Setup

Install by cloning this repo, cd'ing into it and doing `pip install .[dev]`. Run tests by calling
`pytest`.

The database is contained in a directory that must contain `reference` and `sequences`
subdirectories. `sequences` is where sequences that are added to the database get stored. Sequences
and metadata for each isolate gets put in separate fasta and json files. Although this is slightly
cumbersome when it comes to reading in say all the sequences in a database it was done in order to
be able to easily track changes to sequences and metadata for individual isolates using git. In
practice, file i/o rarely takes longer than a second or two even for 10s of thousands of sequences.
`reference` contains fasta files for each segment that each contain just a single sequence:

```
reference/
├── HA.fasta
├── MP.fasta
├── NA.fasta
├── NP.fasta
├── NS.fasta
├── PA.fasta
├── PB1.fasta
└── PB2.fasta
```

Sequences  added to the database are aligned against the appropriate reference sequence and
trimmed.

## Adding sequences

`fs_addseqs` is used to add sequences to the database:

```
usage: fs_addseqs [-h] [--db_dir DB_DIR] new_data_dir

Add sequences to a fluseq db.

positional arguments:
  new_data_dir     Directory containing new data to add. Must contain one or
                   more fasta files and `header_pattern.txt`

options:
  -h, --help       show this help message and exit
  --db_dir DB_DIR  Root directory of a fluseq db. (default='.')
```

`header_pattern.txt` should contain a regex that has named capture groups for extracting metadata
from fasta file headers. E.g. a file for a pipe delimited fasta header might look like this:

```
(?P<isolate_id>.*?)\|
(?P<dna_insdc>.*?)\|
(?P<segment>.*?)\|
(?P<segment_number>.*?)\|
(?P<isolate_name>.*?)\|
(?P<type>.*?)\|
(?P<collection_date>.*?)\|
(?P<identifier>.*?)\|
(?P<dna_accession>.*?)\|
(?P<clade>.*?)\|
(?P<passage>.*?)\|
(?P<lineage>)
```

Sequences that are added get put in `sequences/{segment}/{isolate_id}.fasta`. Metadata gets put in
`sequences/{segment}/{isolate_id}.json`

## Adding additional metadata

`fs_addmeta` can be used to add additional metadata. Warnings are raised if metadata is added for
isolates that do not have sequences in the database. Errors are thrown if data is added that
contradicts existing metadata.

```
usage: fs_addmeta [-h] [--db_dir DB_DIR] table

Add metadata to a fluseq db.

positional arguments:
  table            CSV file containing metadata to add. First column contains
                   isolate_id. Remaining columns are added to the isolate_id's
                   data.

options:
  -h, --help       show this help message and exit
  --db_dir DB_DIR  Root directory of a fluseq db. Default='.'.
```

## Dumping data

`fs_dumpmeta` prints metadata in the database as a CSV:

```
usage: fs_dumpmeta [-h] [--db_dir DB_DIR]

Dump metadata as CSV to stdout

options:
  -h, --help       show this help message and exit
  --db_dir DB_DIR  Root directory of a fluseq db. Default='.'.
```

`fs_dumpseq` prints fasta sequences to stdout:

```
usage: fs_dumpseq [-h] [--id_file ID_FILE] --segment SEGMENT [--db_dir DB_DIR]

Dump sequences from a fluseq db.

options:
  -h, --help         show this help message and exit
  --id_file ID_FILE  Optional file of sequence IDs to print, one per line.
  --segment SEGMENT  Print sequences of this segment.
  --db_dir DB_DIR    Root directory of a fluseq db. Default='.'.
```

`fs_dumpseq_clade` prints fasta sequences of a single clade. The clade is defined by first
rooting the tree using an outgroup and then printing all descendants of the most recent common
ancestor (MRCA) of isolates passed as the ingroup. This program requires building a tree, so can be
slow depending on the size of the database. The tree search uses `fasttree -nt -gtr -nosupport`.

```
usage: fs_dumpseq_clade [-h] [--segments SEGMENTS [SEGMENTS ...]] --outgroup
                        OUTGROUP [--ingroup INGROUP [INGROUP ...]]
                        [--out_dir OUT_DIR] [--db_dir DB_DIR]

Write fasta file comprising a clade.

options:
  -h, --help            show this help message and exit
  --segments SEGMENTS [SEGMENTS ...]
                        Segment(s) to identify clades in. If multiple segments
                        are passed then the union of clade members is output.
                        If this argument is omitted, then all segments are
                        used.
  --outgroup OUTGROUP   Tree is rooted on the outgroup.
  --ingroup INGROUP [INGROUP ...]
                        Descendants of the MRCA of the ingroup are included in
                        the output.
  --out_dir OUT_DIR     Directory to write fasta files for each segment.
  --db_dir DB_DIR       Root directory of a fluseq db. Default='.'.
```
