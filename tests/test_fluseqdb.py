import pytest
import fluseqdb as fsdb


class TestParseHeader:

    def test_a(self):
        pattern = (
            r"(?P<isolate_id>.*?)\|"
            r"(?P<dna_insdc>.*?)\|"
            r"(?P<segment>.*?)\|"
            r"(?P<segment_number>.*?)\|"
            r"(?P<isolate_name>.*?)\|"
            r"(?P<type>.*?)\|"
            r"(?P<collection_date>.*?)\|"
            r"(?P<identifier>.*?)\|"
            r"(?P<dna_accession>.*?)\|"
            r"(?P<clade>.*?)\|"
            r"(?P<passage>.*?)\|"
            r"(?P<lineage>)"
        )
        header = "EPI_ISL_19014384||NP|5|A/dairy_cow/Texas/24-008749-001/2024|A_/_H5N1|2024-03-20|NP_A/dairy_cow/Texas/24-008749-001-original/2024(H5N1)|EPI3158664|2.3.4.4b|Original|"
        metadata = fsdb.add.parse_header(header=header, pattern=pattern)
        assert metadata["segment"] == "NP"

    def test_b(self):
        pattern = (
            r"(?P<isolate_id>.*?)\|"
            r"(?P<dna_insdc>.*?)\|"
            r"(?P<segment>.*?)\|"
            r"(?P<segment_number>.*?)\|"
            r"(?P<isolate_name>.*?)\|"
            r"(?P<type>.*?)\|"
            r"(?P<collection_date>.*?)\|"
            r"(?P<identifier>.*?)\|"
            r"(?P<dna_accession>.*?)\|"
            r"(?P<clade>.*?)\|"
            r"(?P<passage>.*?)\|"
            r"(?P<lineage>)"
        )
        header = "EPI_ISL_19091706||NA|6|A/dairy_cow/New_Mexico/A240920343-108/2024|A_/_H5N1|2024-04-04|A240920343-108|NA|EPI3257806|2.3.4.4b|MDCK1|"
        metadata = fsdb.add.parse_header(header=header, pattern=pattern)
        assert metadata["segment"] == "NA"

    def test_c(self):
        pattern = r"^(?P<isolate_id>.*?)\|(?P<segment>.*?)$"
        header = "MDCKSaleble2_S12|PA"
        metadata = fsdb.add.parse_header(header=header, pattern=pattern)
        assert metadata["segment"] == "PA"


class TestFirstLastNonGap:

    def test_simple(self):
        assert (3, 6) == fsdb.fluseqdb.first_and_last_non_gap("---ACTG---")

    def test_internal_gap(self):
        assert (3, 8) == fsdb.fluseqdb.first_and_last_non_gap("---AC--TG---")


class TestAlignToReference:

    def test_matching_seqs(self):
        output = fsdb.align_to_reference(reference_seq="ACGT", input_seq="ACGT")
        assert output == "ACGT"

    def test_reference_shorter(self):
        output = fsdb.align_to_reference(reference_seq="ACCG", input_seq="TTTACCGTTT")
        assert output == "ACCG"

    def test_input_shorter(self):
        output = fsdb.align_to_reference(reference_seq="ACTG", input_seq="ACG")
        assert output == "AC-G"

    def test_internal_gaps_in_reference_raise_error(self):
        with pytest.raises(ValueError):
            fsdb.align_to_reference(reference_seq="ACTACT", input_seq="ACTGACT")
