"""Golden tests for the accession -> refget contract shared with the MaveDB API.

Why this exists
---------------
A VRS allele digest covers the ``refgetAccession`` of the sequence its location sits on. The mapper and
MaveDB's reverse-translation (RT) job each build alleles, and the allele table deduplicates on that
digest. If the two resolve an accession to different sequences they mint different digests for the same
variant, and the copies never merge.

That happened: the mapper's fetcher chain lacked hgvs's SeqRepo-backed ``SeqFetcher``, so cdot returned a
transcript assembled from the genome (7207 nt) rather than NCBI's NM_007294.3 record (7224 nt), and the
mapper minted ``SQ.bh0R…`` where RT minted ``SQ.jj1R…``. Nothing failed: both digests were internally
consistent, and only reverse-translation's fold-in and the convergent/projection labelling broke.

The two repos run different VRS versions and cannot import each other, so this contract is pinned by a
fixture instead of a shared module. ``canonical_accession_sequences.json`` holds NCBI's records and the
refget each must have. The API repo carries a byte-identical copy and asserts the same thing against its
own implementation. If either repo drifts, its copy of these tests fails.

Expected refgets are ``sha512t24u(sequence)``, not opaque golden strings, so they are correct by
construction and independent of either implementation.
"""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from biocommons.seqrepo import SeqRepo
from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from ga4gh.core import sha512t24u
from hgvs.dataproviders.seqfetcher import SeqFetcher
from hgvs.exceptions import HGVSDataNotAvailableError

from dcd_mapping import lookup, vrs_map
from dcd_mapping.exceptions import (
    AmbiguousReferenceSequenceError,
    ReferenceSequenceNotFoundError,
    ReferenceSequenceProvisioningError,
)

FIXTURE = Path(__file__).parent / "fixtures" / "canonical_accession_sequences.json"
ACCESSIONS = json.loads(FIXTURE.read_text())["accessions"]


@pytest.fixture
def seqrepo_dir(tmp_path: Path) -> Path:
    """Build a small SeqRepo of only the canonical fixture sequences, stored under ``refseq:<accession>``."""
    root = tmp_path / "seqrepo"
    root.mkdir()
    sr = SeqRepo(str(root), writeable=True)
    for accession, record in ACCESSIONS.items():
        sr.store(record["sequence"], [{"namespace": "refseq", "alias": accession}])
    sr.commit()
    return root


@pytest.fixture
def seqrepo_access(seqrepo_dir: Path) -> SeqRepoAccess:
    return SeqRepoAccess(SeqRepo(str(seqrepo_dir)))


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_fixture_refget_is_digest_of_sequence(accession: str):
    """Guards the fixture itself against a hand-edit that desynchronises sequence and refget."""
    record = ACCESSIONS[accession]
    assert record["expected_refget"] == f"SQ.{sha512t24u(record['sequence'].encode())}"
    assert record["expected_refget"] not in record["non_canonical_refgets"]


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_resolve_refget_is_the_canonical_refget(accession, seqrepo_access):
    record = ACCESSIONS[accession]
    refget = lookup.resolve_refget(accession, seqrepo_access)
    assert refget == record["expected_refget"]
    assert refget not in record["non_canonical_refgets"]


def test_missing_accession_fails_closed(seqrepo_access):
    with pytest.raises(ReferenceSequenceNotFoundError):
        lookup.resolve_refget("NM_000000.1", seqrepo_access)


def test_more_than_one_refget_is_an_error_not_a_pick():
    """The upstream alias lookup returns a set, so taking the first of several is order-dependent."""
    seqrepo = MagicMock()
    seqrepo.translate_sequence_identifier.return_value = [
        "ga4gh:SQ.aaa",
        "ga4gh:SQ.bbb",
    ]
    with pytest.raises(AmbiguousReferenceSequenceError):
        lookup.resolve_refget("NM_007294.3", seqrepo)


def test_repeated_alias_for_one_refget_is_not_ambiguous():
    seqrepo = MagicMock()
    seqrepo.translate_sequence_identifier.return_value = [
        "ga4gh:SQ.aaa",
        "ga4gh:SQ.aaa",
    ]
    assert lookup.resolve_refget("NM_007294.3", seqrepo) == "SQ.aaa"


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_mapping_an_accession_neither_writes_to_seqrepo_nor_asks_cdot(
    accession, seqrepo_access, seqrepo_dir
):
    """Regression for the original defect: a mapping run used to fetch the sequence from cdot and store it
    under the accession alias. It must now only read, so cdot's sequence can never displace NCBI's.
    """
    record = ACCESSIONS[accession]
    cdot = MagicMock()
    cdot.get_seq.return_value = (
        "ACGT" * 100
    )  # any sequence that is not the canonical record
    with (
        patch.object(lookup, "get_seqrepo", return_value=seqrepo_access),
        patch.object(lookup, "cdot_rest", return_value=cdot),
    ):
        vrs_map.ensure_accession_in_seqrepo(accession)

    cdot.get_seq.assert_not_called()
    reopened = SeqRepoAccess(SeqRepo(str(seqrepo_dir)))
    assert lookup.resolve_refget(accession, reopened) == record["expected_refget"]


def test_mapping_an_accession_absent_from_seqrepo_fails_closed(seqrepo_access):
    with (
        patch.object(lookup, "get_seqrepo", return_value=seqrepo_access),
        pytest.raises(ReferenceSequenceNotFoundError),
    ):
        vrs_map.ensure_accession_in_seqrepo("NM_000000.1")


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_seqfetcher_serves_ncbi_records_ahead_of_the_fasta_fallback(
    accession, seqrepo_dir, monkeypatch
):
    """The fetcher chain feeding hgvs validation and normalisation must resolve accessions from SeqRepo
    before the genome-assembled FASTA fallback. The FASTA fetcher is stubbed to refuse, so a chain
    without the SeqRepo-backed fetcher in front (the original defect) cannot serve the sequence.
    """
    monkeypatch.setenv("HGVS_SEQREPO_DIR", str(seqrepo_dir))
    fasta = MagicMock(source="stub fasta")
    fasta.fetch_seq.side_effect = HGVSDataNotAvailableError("no fasta in test")
    with (
        patch.object(lookup, "FastaSeqFetcher", return_value=fasta),
        patch.object(lookup, "GENOMIC_FASTA_FILES", ["stub.fna"]),
    ):
        chain = lookup.seqfetcher()

    assert isinstance(chain.seq_fetchers[0], SeqFetcher)
    sequence = chain.fetch_seq(accession)
    record = ACCESSIONS[accession]
    assert sequence == record["sequence"]
    assert f"SQ.{sha512t24u(sequence.encode())}" == record["expected_refget"]


# --- Shared-store writes: short transactions and Ensembl provisioning -----------------------------


def _shared_writer(seqrepo_dir: Path) -> SeqRepoAccess:
    """Open the SeqRepo the way the mapper does: writeable, with the shared-write busy timeout."""
    sr = SeqRepo(str(seqrepo_dir), writeable=True)
    lookup.configure_seqrepo_for_shared_writes(sr)
    return SeqRepoAccess(sr)


def _reader(seqrepo_dir: Path) -> SeqRepoAccess:
    """Open a fresh read-only connection, as reverse translation does."""
    return SeqRepoAccess(SeqRepo(str(seqrepo_dir)))


ENST = "ENST00000999999.1"
TRANSCRIPT = "ACGT" * 60


def _refget(sequence: str) -> str:
    return f"SQ.{sha512t24u(sequence.encode())}"


def test_busy_timeout_is_raised_for_shared_writes(seqrepo_dir):
    writer = _shared_writer(seqrepo_dir)
    for db in (writer.sr.aliases._db, writer.sr.sequences._db):
        assert (
            db.execute("PRAGMA busy_timeout").fetchone()[0]
            == lookup.SEQREPO_BUSY_TIMEOUT_MS
        )


def test_store_sequence_commits_so_it_does_not_hold_the_write_lock(seqrepo_dir):
    """A mapping process that stored a target sequence used to keep an open transaction, blocking every
    other writer for the life of the process. A second writer must now get through at once.
    """
    mapper = _shared_writer(seqrepo_dir)
    with patch.object(vrs_map, "get_seqrepo", return_value=mapper):
        sequence_id = vrs_map.store_sequence("ACGT" * 50)

    other = SeqRepo(str(seqrepo_dir), writeable=True)
    other.aliases._db.execute("PRAGMA busy_timeout = 200")
    other.sequences._db.execute("PRAGMA busy_timeout = 200")
    other.store("TTTT" * 50, [{"namespace": "test", "alias": "second-writer"}])
    other.commit()

    assert list(
        _reader(seqrepo_dir).sr.aliases.find_aliases(
            namespace="ga4gh", alias=sequence_id
        )
    )


def test_provisioning_an_ensembl_transcript_persists_it_for_other_readers(seqrepo_dir):
    writer = _shared_writer(seqrepo_dir)
    builds = {"GRCh38.fna": TRANSCRIPT, "GRCh37.fna": TRANSCRIPT}
    with (
        patch.object(lookup, "get_seqrepo", return_value=writer),
        patch.object(lookup, "_assemble_ensembl_transcript", return_value=builds),
    ):
        refget = lookup.provision_ensembl_transcript(ENST)

    assert refget == _refget(TRANSCRIPT)
    assert lookup.resolve_refget(ENST, _reader(seqrepo_dir)) == refget


def test_provisioning_refuses_when_genome_builds_disagree(seqrepo_dir):
    writer = _shared_writer(seqrepo_dir)
    builds = {"GRCh38.fna": TRANSCRIPT, "GRCh37.fna": TRANSCRIPT[:-1] + "A"}
    with (
        patch.object(lookup, "get_seqrepo", return_value=writer),
        patch.object(lookup, "_assemble_ensembl_transcript", return_value=builds),
        pytest.raises(ReferenceSequenceProvisioningError),
    ):
        lookup.provision_ensembl_transcript(ENST)
    with pytest.raises(ReferenceSequenceNotFoundError):
        lookup.resolve_refget(ENST, _reader(seqrepo_dir))


@pytest.mark.parametrize(
    "builds",
    [{}, {"GRCh38.fna": ""}, {"GRCh38.fna": "ACGTXACGT"}],
    ids=["no-build", "empty", "not-nucleotides"],
)
def test_provisioning_refuses_an_unusable_assembly(seqrepo_dir, builds):
    writer = _shared_writer(seqrepo_dir)
    with (
        patch.object(lookup, "get_seqrepo", return_value=writer),
        patch.object(lookup, "_assemble_ensembl_transcript", return_value=builds),
        pytest.raises(ReferenceSequenceProvisioningError),
    ):
        lookup.provision_ensembl_transcript(ENST)


def test_a_missing_ensembl_transcript_is_provisioned_when_mapping(seqrepo_dir):
    writer = _shared_writer(seqrepo_dir)
    builds = {"GRCh38.fna": TRANSCRIPT}
    with (
        patch.object(lookup, "get_seqrepo", return_value=writer),
        patch.object(lookup, "_assemble_ensembl_transcript", return_value=builds),
    ):
        vrs_map.ensure_accession_in_seqrepo(ENST)
    assert lookup.resolve_refget(ENST, _reader(seqrepo_dir)) == _refget(TRANSCRIPT)


def test_a_present_ensembl_transcript_is_never_reassigned(seqrepo_dir):
    writer = _shared_writer(seqrepo_dir)
    writer.sr.store(TRANSCRIPT, [{"namespace": "ensembl", "alias": ENST}])
    writer.sr.commit()
    assemble = MagicMock(return_value={"GRCh38.fna": "TTTT" * 60})
    with (
        patch.object(lookup, "get_seqrepo", return_value=writer),
        patch.object(lookup, "_assemble_ensembl_transcript", assemble),
    ):
        vrs_map.ensure_accession_in_seqrepo(ENST)

    assemble.assert_not_called()
    assert lookup.resolve_refget(ENST, _reader(seqrepo_dir)) == _refget(TRANSCRIPT)


def test_a_missing_refseq_transcript_is_never_assembled(seqrepo_dir):
    """Only Ensembl transcripts are exact when assembled from the genome; RefSeq records can differ."""
    writer = _shared_writer(seqrepo_dir)
    assemble = MagicMock()
    with (
        patch.object(lookup, "get_seqrepo", return_value=writer),
        patch.object(lookup, "_assemble_ensembl_transcript", assemble),
        pytest.raises(ReferenceSequenceNotFoundError),
    ):
        vrs_map.ensure_accession_in_seqrepo("NM_000000.1")
    assemble.assert_not_called()
