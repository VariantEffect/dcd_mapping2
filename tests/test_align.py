"""Test ``align`` module.

Todo:
----
* Mock the BLAT call/result file


"""

from unittest.mock import patch

import numpy as np
import pytest
from Bio import Align as BioAlign
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dcd_mapping.align import (
    _blat_score,
    _build_alignment_qc,
    _compact_alignment_string,
    _get_best_hsp,
    align,
)
from dcd_mapping.exceptions import AlignmentError
from dcd_mapping.schemas import (
    AlignmentQc,
    AlignmentResult,
    ScoresetMetadata,
    TargetGene,
    TargetType,
)


def _make_aln(
    query_seq: str,
    target_seq: str,
    q_start: int = 0,
    t_start: int = 0,
    coords: "np.ndarray | None" = None,
) -> BioAlign.Alignment:
    """Build a minimal Bio.Align.Alignment for unit-testing _build_alignment_qc."""
    n = len(query_seq)
    t_len = max(t_start + n + 10, len(target_seq) + 10)
    # Pad target so sequence slicing is always in bounds
    padded_target = target_seq.ljust(t_len, "N")
    target = SeqRecord(Seq(padded_target), id="chr1")
    query = SeqRecord(Seq(query_seq), id="query")
    if coords is None:
        coords = np.array([[t_start, t_start + n], [q_start, q_start + n]])
    return BioAlign.Alignment(sequences=[target, query], coordinates=coords)


class TestBuildAlignmentQc:
    def test_exact_match(self):
        seq = "ACGTACGTAC"
        aln = _make_aln(seq, seq)
        qc = _build_alignment_qc(aln)
        assert isinstance(qc, AlignmentQc)
        assert qc.mismatch_count == 0
        assert qc.gap_count == 0
        assert qc.mismatch_positions == []
        assert qc.mismatch_positions_unavailable is False
        assert qc.gap_intervals == []
        assert qc.cigar == "10M"
        assert qc.percent_identity == pytest.approx(100.0)

    def test_single_mismatch(self):
        query = "ACGTACGTAC"
        # position 3: T->G mismatch
        target = "ACGGACGTAC"
        aln = _make_aln(query, target)
        qc = _build_alignment_qc(aln)
        assert qc.mismatch_count == 1
        assert qc.mismatch_positions == [3]

    def test_gap_in_ref(self):
        """A gap appears when target has fewer bases than query (insertion in query)."""
        target_seq = "ACGTAACGTA"  # 10 bases
        query_seq = "ACGTAGGACGTA"  # 12 bases (2 extra at position 5-7)
        target = SeqRecord(Seq(target_seq + "N" * 10), id="chr1")
        query = SeqRecord(Seq(query_seq), id="query")
        coords = np.array([[0, 5, 5, 10], [0, 5, 7, 12]])
        aln = BioAlign.Alignment(sequences=[target, query], coordinates=coords)
        qc = _build_alignment_qc(aln)
        assert qc.gap_count == 1
        assert len(qc.gap_intervals) == 1

    def test_gap_in_query(self):
        """A gap appears when query has fewer bases than target (deletion in query)."""
        target_seq = "ACGTANNNACGTA"  # 13 bases
        query_seq = "ACGTAACGTA"  # 10 bases
        target = SeqRecord(Seq(target_seq + "N" * 10), id="chr1")
        query = SeqRecord(Seq(query_seq), id="query")
        coords = np.array([[0, 5, 8, 13], [0, 5, 5, 10]])
        aln = BioAlign.Alignment(sequences=[target, query], coordinates=coords)
        qc = _build_alignment_qc(aln)
        assert qc.gap_count == 1
        assert len(qc.gap_intervals) == 1
        gs, ge = qc.gap_intervals[0]
        assert ge - gs == 3

    def test_multi_block_alignment(self):
        """Two aligned blocks with a gap between them."""
        seq = "ACGTAACGTA"  # 10-base query
        target_seq = seq[:5] + "NNN" + seq[5:]  # 13-base target with 3-base insertion
        target = SeqRecord(Seq(target_seq + "N" * 10), id="chr1")
        query = SeqRecord(Seq(seq), id="query")
        coords = np.array([[0, 5, 8, 13], [0, 5, 5, 10]])
        aln = BioAlign.Alignment(sequences=[target, query], coordinates=coords)
        qc = _build_alignment_qc(aln)
        # alignment_length = target span = coords[0][-1] - coords[0][0] = 13
        assert qc.alignment_length == 13
        assert qc.gap_count == 1
        assert qc.mismatch_count == 0

    def test_undefined_sequence_fallback(self):
        """_build_alignment_qc falls back gracefully when sequence content is unavailable.

        Biopython's PSL parser leaves SeqRecords with undefined content when BLAT's
        pslx output omits per-block sequences (common in protein-vs-protein mode).
        Slicing raises UndefinedSequenceError; the function should fall back to an
        empty mismatch_positions list rather than crashing, and must set
        mismatch_positions_unavailable=True so callers know at_mismatched_locus is
        unreliable when mismatch_count > 0.
        """
        target = SeqRecord(Seq(None, length=20), id="reference")
        query = SeqRecord(Seq(None, length=10), id="query")
        coords = np.array([[0, 10], [0, 10]])
        aln = BioAlign.Alignment(sequences=[target, query], coordinates=coords)

        qc = _build_alignment_qc(aln)

        assert qc.mismatch_positions == []
        assert qc.mismatch_positions_unavailable is True
        # gap_intervals and gap_count are coordinate-driven (no sequence needed)
        assert qc.gap_count == 0

    def test_mismatch_positions_after_sequence_attachment(self):
        """Attaching sequences to an undefined alignment enables per-base mismatch detection.

        This mirrors what align_target_to_protein does before calling _build_alignment_qc:
        the PSL gives coordinates; we supply the original sequences so that slicing works.
        """
        reference = "MKALVDQSLR"
        query_seq = "MKALXDQZLR"  # mismatches at positions 4 ('X') and 7 ('Z')

        target = SeqRecord(Seq(None, length=len(reference)), id="reference")
        query = SeqRecord(Seq(None, length=len(query_seq)), id="query")
        coords = np.array([[0, len(reference)], [0, len(query_seq)]])
        aln = BioAlign.Alignment(sequences=[target, query], coordinates=coords)

        # Mimic align_target_to_protein: attach sequences before QC
        aln.target.seq = Seq(reference)
        aln.query.seq = Seq(query_seq)

        qc = _build_alignment_qc(aln)

        assert 4 in qc.mismatch_positions
        assert 7 in qc.mismatch_positions
        assert len(qc.mismatch_positions) == 2
        assert qc.mismatch_count == 2

    def test_minus_strand_alignment(self):
        """_build_alignment_qc works for minus-strand PSL alignments.

        Biopython's SAM formatter raises 'Unequal step sizes' for PSL
        alignments where the query is on the minus strand, because query
        coordinates are stored in descending order.  We build CIGAR from
        coordinates directly so this is a non-issue.
        """
        # 8-base sequences; minus-strand: query[8:0] aligns to target[10:18].
        # Biopython stores minus-strand query coords as descending.
        seq = "ACGTACGT"
        target = SeqRecord(Seq(seq * 20), id="chr1")
        query = SeqRecord(Seq(seq), id="query")
        coords = np.array([[10, 18], [8, 0]])
        aln = BioAlign.Alignment(sequences=[target, query], coordinates=coords)

        qc = _build_alignment_qc(aln)

        assert qc.cigar == "8M"
        assert qc.gap_count == 0


class TestCompactAlignmentString:
    """Unit tests for _compact_alignment_string.

    The function collapses consecutive all-'?' target-sequence blocks (intronic
    regions in cDNA→genome BLAT output) into summary lines of the form
    ``... [N bp gap: chrom:start-end] ...``.
    """

    def _make_group(
        self, chrom: str, start: int, seq: str, match_row: str, query_seq: str
    ) -> str:
        """Build a single Biopython-style 3-line display group."""
        return f"{chrom} {start} {seq}\n           {match_row}\n query 0    {query_seq}"

    def test_noop_on_no_gap_blocks(self):
        """A string with no all-'?' target rows is returned unchanged."""
        group = self._make_group("chr1", 100, "ACGT", "||||", "ACGT")
        assert _compact_alignment_string(group) == group

    def test_single_gap_block_collapsed(self):
        """One all-'?' block is replaced with a summary line."""
        exon1 = self._make_group("chr1", 100, "ACGT", "||||", "ACGT")
        intron = self._make_group("chr1", 104, "????", "    ", "----")
        exon2 = self._make_group("chr1", 108, "TTTT", "||||", "TTTT")
        inp = "\n\n".join([exon1, intron, exon2])
        out = _compact_alignment_string(inp)
        # Intron block replaced; exon blocks preserved
        assert "... [4 bp gap: chr1:104-108] ..." in out
        assert "ACGT" in out
        assert "TTTT" in out
        # No raw '????' remaining
        assert "????" not in out

    def test_consecutive_gap_blocks_merged(self):
        """Multiple consecutive all-'?' blocks are merged into a single summary."""
        exon = self._make_group("chr1", 0, "ACGT", "||||", "ACGT")
        intron1 = self._make_group("chr1", 4, "?" * 60, " " * 60, "-" * 60)
        intron2 = self._make_group("chr1", 64, "?" * 60, " " * 60, "-" * 60)
        inp = "\n\n".join([exon, intron1, intron2])
        out = _compact_alignment_string(inp)
        # Both intron blocks merged into one summary
        assert "... [120 bp gap: chr1:4-124] ..." in out
        assert out.count("... [") == 1

    def test_mixed_sequence_block_not_collapsed(self):
        """A block with mixed '?' and real bases is left as-is."""
        mixed = self._make_group("chr1", 0, "AC?T", "|| |", "ACGT")
        out = _compact_alignment_string(mixed)
        assert "AC?T" in out
        assert "... [" not in out

    def test_protein_alignment_noop(self):
        """Protein alignments have no '?' characters; function is a no-op."""
        protein_group = "NP_001234.1 0 MAAAA\n          |||||\n query 0   MAAAA"
        assert _compact_alignment_string(protein_group) == protein_group

    def test_comma_formatting_in_large_gap(self):
        """Gap length ≥ 1000 uses comma-formatted number."""
        exon = self._make_group("chr1", 0, "ACGT", "||||", "ACGT")
        # Build consecutive 60-bp intron blocks totalling > 1000 bp
        introns = [
            self._make_group("chr1", 4 + i * 60, "?" * 60, " " * 60, "-" * 60)
            for i in range(17)  # 17 * 60 = 1020 bp
        ]
        inp = "\n\n".join([exon, *introns])
        out = _compact_alignment_string(inp)
        assert "1,020 bp gap" in out


class TestBlatScore:
    """_blat_score and _get_best_hsp regression tests.

    The key invariant: _get_best_hsp must use BLAT PSL score
    (identities - mismatches - qNumInsert - tNumInsert), NOT raw identity
    count.  These tests construct two alignments where the two metrics
    disagree and verify that BLAT-score-based selection wins.
    """

    def _make_noisy_aln(self) -> BioAlign.Alignment:
        """Alignment A: 10 identity matches + 8 mismatches in one block.

        BLAT score  = 10 - 8 = 2
        identity count = 10
        """
        # First 10 chars identical, last 8 mismatched (N vs T)
        target_seq = "ACGTACGTAC" + "TTTTTTTT"
        query_seq = "ACGTACGTAC" + "NNNNNNNN"
        padded = target_seq.ljust(len(target_seq) + 10, "A")
        target = SeqRecord(Seq(padded), id="chr1")
        query = SeqRecord(Seq(query_seq), id="query_noisy")
        coords = np.array([[0, len(query_seq)], [0, len(query_seq)]])
        return BioAlign.Alignment(sequences=[target, query], coordinates=coords)

    def _make_clean_aln(self) -> BioAlign.Alignment:
        """Alignment B: 7 perfect matches, no mismatches.

        BLAT score  = 7 - 0 = 7
        identity count = 7

        Shorter but cleaner: should be preferred by _get_best_hsp.
        """
        seq = "ACGTACG"
        padded = seq.ljust(len(seq) + 10, "A")
        target = SeqRecord(Seq(padded), id="chr1")
        query = SeqRecord(Seq(seq), id="query_clean")
        coords = np.array([[0, len(seq)], [0, len(seq)]])
        return BioAlign.Alignment(sequences=[target, query], coordinates=coords)

    def test_blat_score_noisy_lower_than_clean(self):
        """Noisy alignment has more identities but a lower BLAT score."""
        noisy = self._make_noisy_aln()
        clean = self._make_clean_aln()

        noisy_counts = noisy.counts()
        clean_counts = clean.counts()

        # Identity count: noisy wins
        assert noisy_counts.identities > clean_counts.identities, (
            "Noisy alignment must have more raw identities for the test to be meaningful"
        )
        # BLAT score: clean wins
        assert _blat_score(noisy) < _blat_score(clean), (
            "Clean alignment must have the higher BLAT score for the test to be meaningful"
        )

    def test_get_best_hsp_prefers_blat_score_over_identity_count(self):
        """_get_best_hsp selects the alignment with the higher BLAT score, not
        the higher raw identity count.

        Regression guard: before the fix, the selector used
        ``a.counts().identities``, which would have returned the noisy
        alignment.  After the fix it must return the clean one.
        """
        noisy = self._make_noisy_aln()
        clean = self._make_clean_aln()

        best = _get_best_hsp([noisy, clean])

        assert best is clean, (
            "_get_best_hsp chose the noisy alignment (higher identity count) "
            "instead of the clean one (higher BLAT score). "
            "The ranking criterion has regressed to identity-count-only."
        )


def check_alignment_result_equality(actual: AlignmentResult, expected: AlignmentResult):
    assert actual.chrom == expected.chrom
    assert actual.strand == expected.strand
    assert actual.coverage == pytest.approx(expected.coverage)
    assert actual.percent_identity == pytest.approx(expected.percent_identity)
    assert actual.query_range.start == expected.query_range.start
    assert actual.query_range.end == expected.query_range.end
    for a, e in zip(actual.query_subranges, expected.query_subranges, strict=False):
        assert a.start == e.start
        assert a.end == e.end
    assert len(actual.query_subranges) == len(expected.query_subranges)
    assert actual.hit_range.start == expected.hit_range.start
    assert actual.hit_range.end == expected.hit_range.end
    for a, e in zip(actual.hit_subranges, expected.hit_subranges, strict=False):
        assert a.start == e.start
        assert a.end == e.end
    assert len(actual.hit_subranges) == len(expected.hit_subranges)


def test_align_src_catalytic_domain(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000041-a-1"""
    urn = "urn:mavedb:00000041-a-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = AlignmentResult(**align_result_fixture[urn])
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_hbb(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000018-a-1"""
    urn = "urn:mavedb:00000018-a-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = AlignmentResult(**align_result_fixture[urn])
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_ube2i(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` on urn:mavedb:00000001-a-4"""
    urn = "urn:mavedb:00000001-a-4"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = AlignmentResult(**align_result_fixture[urn])
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_scn5a(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000098-a-1"""
    urn = "urn:mavedb:00000098-a-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = AlignmentResult(**align_result_fixture[urn])
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_raf(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000061-h-1."""
    urn = "urn:mavedb:00000061-h-1"
    scoreset_metadata = scoreset_metadata_fixture[urn]
    align_result = align(scoreset_metadata)
    expected = AlignmentResult(**align_result_fixture[urn])
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_tp53(scoreset_metadata_fixture, align_result_fixture):
    """Test ``align()`` method on urn:mavedb:00000068-a-1"""
    urn = "urn:mavedb:00000068-a-1"
    metadata = scoreset_metadata_fixture[urn]
    align_result = align(metadata)
    expected = AlignmentResult(**align_result_fixture[urn])
    assert align_result
    check_alignment_result_equality(align_result, expected)


def test_align_raises_on_ambiguous_blat_id():
    """align() must raise AlignmentError when a BLAT result ID matches more than one
    target gene name after stripping whitespace and punctuation.

    Regression guard for the silent-pass bug where ``matches > 1`` constructed
    an error message but never raised it.
    """
    # Two gene names that both normalise to the same BLAT identifier "GENE".
    # "GENE (isoform1)" -> split()[0].replace("(","").replace(")","") -> "GENE"
    # "GENE (isoform2)" -> same result -> ambiguous.
    metadata = ScoresetMetadata(
        urn="urn:mavedb:00000001-a-1",
        target_genes={
            "GENE (isoform1)": TargetGene(
                target_gene_name="GENE (isoform1)",
                target_gene_category=TargetType.PROTEIN_CODING,
                target_sequence="ACGT",
                target_sequence_type="dna",
                target_accession_id=None,
            ),
            "GENE (isoform2)": TargetGene(
                target_gene_name="GENE (isoform2)",
                target_gene_category=TargetType.PROTEIN_CODING,
                target_sequence="ACGT",
                target_sequence_type="dna",
                target_accession_id=None,
            ),
        },
    )
    # _get_blat_output returns a grouped dict keyed by the shortened BLAT id.
    # We fake a single hit for the shared prefix so the disambiguation loop runs.
    fake_grouped = {"GENE": []}
    fake_params = {"aligner": "blat"}

    with (
        patch(
            "dcd_mapping.align._get_blat_output",
            return_value=(fake_grouped, fake_params),
        ),
        pytest.raises(AlignmentError, match="matches multiple target gene names"),
    ):
        align(metadata)
