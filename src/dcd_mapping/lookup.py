"""Handle API lookups to external (non-MaveDB) services.

Data sources/handlers include:

* `CoolSeqTool <https://github.com/GenomicMedLab/cool-seq-tool/>`_
* `Gene Normalizer <https://github.com/cancervariants/gene-normalization>`_
* the `VRS-Python Translator tool <https://github.com/ga4gh/vrs-python>`_
* the UniProt web API
"""

import logging
import os
from pathlib import Path
from typing import Any

import hgvs
import httpx
import polars as pl
from biocommons.seqrepo import SeqRepo
from biocommons.seqrepo.seqaliasdb.seqaliasdb import sqlite3
from cdot.hgvs.dataproviders import ChainedSeqFetcher, FastaSeqFetcher, RESTDataProvider
from cool_seq_tool.app import (
    LRG_REFSEQGENE_PATH,
    MANE_SUMMARY_PATH,
    SEQREPO_ROOT_DIR,
    TRANSCRIPT_MAPPINGS_PATH,
    UTA_DB_URL,
    CoolSeqTool,
)
from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.mappers import (
    AlignmentMapper,
    ExonGenomicCoordsMapper,
    ManeTranscript,
)
from cool_seq_tool.schemas import TranscriptPriority
from cool_seq_tool.sources.mane_transcript_mappings import ManeTranscriptMappings
from cool_seq_tool.sources.transcript_mappings import TranscriptMappings
from cool_seq_tool.sources.uta_database import UtaDatabase
from ga4gh.core._internal.models import Extension, Gene
from ga4gh.vrs._internal.models import (
    Allele,
    LiteralSequenceExpression,
    ReferenceLengthExpression,
    SequenceLocation,
    SequenceReference,
)
from ga4gh.vrs.dataproxy import SeqRepoDataProxy, coerce_namespace
from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.utils.hgvs_tools import HgvsTools
from gene.database import create_db
from gene.query import QueryHandler
from gene.schemas import MatchType, SourceName

from dcd_mapping.exceptions import DataLookupError
from dcd_mapping.resource_utils import CDOT_URL, ENSEMBL_API_URL, request_with_backoff
from dcd_mapping.schemas import (
    GeneLocation,
    ManeDescription,
    TargetGene,
)

__all__ = [
    "CoolSeqToolBuilder",
    "GeneNormalizerBuilder",
    "build_ref_identical_allele",
    "get_chromosome_identifier",
    "get_chromosome_identifier_from_vrs_id",
    "get_gene_location",
    "get_gene_symbol",
    "get_mane_transcripts",
    "get_protein_accession",
    "get_seqrepo",
    "get_sequence",
    "get_transcripts",
    "get_ucsc_chromosome_name",
    "get_uniprot_sequence",
    "translate_hgvs_to_vrs",
    "translate_ref_identical_to_vrs",
]
_logger = logging.getLogger(__name__)

# ---------------------------------- Cdot ---------------------------------- #


GENOMIC_FASTA_FILES = [
    "/home/.local/share/dcd_mapping/GCF_000001405.39_GRCh38.p13_genomic.fna.gz",
    "/home/.local/share/dcd_mapping/GCF_000001405.25_GRCh37.p13_genomic.fna.gz",
]


def seqfetcher() -> ChainedSeqFetcher:
    return ChainedSeqFetcher(*[FastaSeqFetcher(file) for file in GENOMIC_FASTA_FILES])


def cdot_rest() -> RESTDataProvider:
    return RESTDataProvider(url=CDOT_URL, seqfetcher=seqfetcher())


# ---------------------------------- Global ---------------------------------- #


class _TimeoutUtaDatabase(UtaDatabase):
    """UtaDatabase subclass with an increased command timeout.

    The upstream default is 60s, which can be insufficient for range queries
    against the public UTA server's ``tx_exon_aln_v`` view.
    """

    async def create_pool(self) -> None:
        """Create connection pool with a 5-minute command timeout."""
        if not self._connection_pool:
            import asyncpg

            self.args = self._get_conn_args()
            try:
                self._connection_pool = await asyncpg.create_pool(
                    min_size=1,
                    max_size=10,
                    max_inactive_connection_lifetime=3,
                    command_timeout=300,
                    host=self.args["host"],
                    port=self.args["port"],
                    user=self.args["user"],
                    password=self.args["password"],
                    database=self.args["database"],
                )
            except asyncpg.InterfaceError as e:
                _logger.error(
                    "While creating connection pool, encountered exception %s", e
                )
                msg = "Could not create connection pool"
                raise Exception(msg) from e


class CoolSeqToolBuilder:
    """Singleton constructor for ``cool-seq-tool`` instance."""

    def __new__(cls) -> CoolSeqTool:
        """Provide ``CoolSeqTool`` instance. Construct it if unavailable.

        This class temporarily includes some very obnoxious reimplementations of
        CoolSeqTool classes due to some changes introduced in VRS-Python 2a6. We should
        try to clean them up.

        :return: singleton instance of CoolSeqTool
        """

        class _AugmentedSeqRepoAccess(SeqRepoAccess):
            def derive_refget_accession(self, ac: str) -> str | None:
                if ac is None:
                    return None

                if ":" not in ac[1:]:
                    # always coerce the namespace if none provided
                    ac = coerce_namespace(ac)

                refget_accession = None
                try:
                    aliases = self.translate_sequence_identifier(ac, namespace="ga4gh")
                except KeyError:
                    _logger.exception("KeyError when getting refget accession: %s", ac)
                else:
                    if aliases:
                        refget_accession = aliases[0].split("ga4gh:")[-1]

                return refget_accession

        class _AugmentedCoolSeqTool(CoolSeqTool):
            def __init__(
                self,
                transcript_file_path: Path = TRANSCRIPT_MAPPINGS_PATH,
                lrg_refseqgene_path: Path = LRG_REFSEQGENE_PATH,
                mane_data_path: Path = MANE_SUMMARY_PATH,
                db_url: str = UTA_DB_URL,
                sr: SeqRepo | None = None,
            ) -> None:
                if not sr:
                    sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
                self.seqrepo_access = _AugmentedSeqRepoAccess(sr)
                self.transcript_mappings = TranscriptMappings(
                    transcript_file_path=transcript_file_path,
                    lrg_refseqgene_path=lrg_refseqgene_path,
                )
                self.mane_transcript_mappings = ManeTranscriptMappings(
                    mane_data_path=mane_data_path
                )
                self.uta_db = _TimeoutUtaDatabase(db_url=db_url)
                self.alignment_mapper = AlignmentMapper(
                    self.seqrepo_access, self.transcript_mappings, self.uta_db
                )
                self.mane_transcript = ManeTranscript(
                    self.seqrepo_access,
                    self.transcript_mappings,
                    self.mane_transcript_mappings,
                    self.uta_db,
                )
                self.ex_g_coords_mapper = ExonGenomicCoordsMapper(
                    self.seqrepo_access,
                    self.uta_db,
                    self.mane_transcript,
                    self.mane_transcript_mappings,
                )

        if not hasattr(cls, "instance"):
            root_dir = os.environ.get(
                "SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest"
            )
            sr = SeqRepo(root_dir, writeable=True)
            cls.instance = _AugmentedCoolSeqTool(sr=sr)

        return cls.instance


def get_seqrepo() -> SeqRepoAccess:
    """Retrieve SeqRepo access instance."""
    cst = CoolSeqToolBuilder()
    return cst.seqrepo_access


class GeneNormalizerBuilder:
    """Singleton constructor for Gene Normalizer instance."""

    def __new__(cls) -> QueryHandler:
        """Provide Gene Normalizer instance. If an instance has already been
        constructed, close its connection and provide a new one.

        :return: singleton instance of ``QueryHandler`` for Gene Normalizer
        """
        if hasattr(cls, "instance"):
            cls.instance.db.close_connection()
            cls.instance = None

        db = create_db()
        cls.instance = QueryHandler(db)
        return cls.instance


def init_hgvs_tools(self, data_proxy=None):  # noqa: ANN202, ANN001
    """Initialize HgvsTools with cdot as data provider"""
    self.parser = hgvs.parser.Parser()
    self.data_proxy = data_proxy
    cdot_provider = cdot_rest()
    self.normalizer = hgvs.normalizer.Normalizer(cdot_provider, validate=True)
    self.variant_mapper = hgvs.variantmapper.VariantMapper(cdot_provider)


class TranslatorBuilder:
    """Singleton constructor for VRS Translator instance."""

    def __new__(cls, data_proxy: SeqRepoDataProxy) -> AlleleTranslator:
        """Provide translator instance. Constructs it if unavailable. Use a new
        ``data_proxy`` instance that contains a given score row's sequence/ID.

        :return: singleton instance of ``AlleleTranslator``
        """
        if not hasattr(cls, "instance"):
            # monkey patch to use cdot instead of UTA as HgvsTools data provider
            HgvsTools.__init__ = init_hgvs_tools
            tr = AlleleTranslator(data_proxy)
            cls.instance = tr
        else:
            cls.instance.data_proxy = data_proxy
        return cls.instance


# ----------------------------------- UTA ----------------------------------- #


async def check_uta() -> None:
    """Check that UTA connection appears to be working.

    :raise LookupError: if schema check fails. Realistically, if UTA isn't working right,
        another error will probably get raised first.
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"select * from {uta.schema}.meta"  # noqa: S608
    result = await uta.execute_query(query)
    if not result:
        msg = "UTA schema check failed. No results returned."
        _logger.error(msg)
        raise DataLookupError(msg)


async def get_protein_accession(transcript: str) -> str | None:
    """Retrieve protein accession for a transcript.

    :param transcript: transcript accession, e.g. ``"NM_002529.3"``
    :return: protein accession if successful
    """
    try:
        uta = CoolSeqToolBuilder().uta_db
        query = f"""
        SELECT pro_ac FROM {uta.schema}.associated_accessions
        WHERE tx_ac = '{transcript}'
        """  # noqa: S608
        result = await uta.execute_query(query)
    except Exception as e:
        _logger.exception(
            "Failed to get protein accession for transcript %s", transcript
        )
        raise DataLookupError from e
    if result:
        return result[0]["pro_ac"]
    return None


async def get_transcripts(
    chromosome_ac: str, start: int, end: int
) -> list[tuple[str, str]]:
    """Get transcript accessions matching given parameters (excluding non-coding RNA),
    returning both the transcript accession and HGNC symbol.

    :param chromosome: chromosome accession (e.g. ``"NC_000007.13"``)
    :param start: starting position
    :param end: ending position
    :return: candidate transcript accessions and HGNC symbols
    """
    try:
        uta = CoolSeqToolBuilder().uta_db
        query = f"""
        SELECT tx_ac, hgnc
        FROM {uta.schema}.tx_exon_aln_v
        WHERE ({start} BETWEEN alt_start_i AND alt_end_i OR {end} BETWEEN alt_start_i AND alt_end_i)
        AND alt_ac = '{chromosome_ac}'
        AND tx_ac NOT LIKE 'NR_%';
        """  # noqa: S608
        result = await uta.execute_query(query)
    except Exception as e:
        _logger.exception(
            "Failed to get transcripts for %s:%d-%d", chromosome_ac, start, end
        )
        raise DataLookupError from e

    return [(row["tx_ac"], row["hgnc"]) for row in result]


# ------------------------------ Gene Normalizer ------------------------------ #


def check_gene_normalizer() -> None:
    q = GeneNormalizerBuilder()
    if (not q.db.check_schema_initialized()) or not (q.db.check_tables_populated()):
        msg = "Gene Normalizer database schema check failed. No results returned."
        _logger.error(msg)
        raise DataLookupError(msg)
    if q.normalize("BRAF").match_type == MatchType.NO_MATCH:
        msg = "Gene Normalizer returned no normalization results for BRAF. This indicates an underlying issue with the database that should be investigated."
        _logger.error(msg)
        raise DataLookupError(msg)


def _get_hgnc_symbol(term: str) -> str | None:
    """Fetch HGNC symbol from gene term.

    :param term: gene referent
    :return: gene symbol if available
    """
    q = GeneNormalizerBuilder()
    result = q.normalize_unmerged(term)
    hgnc = result.source_matches.get(SourceName.HGNC)
    if hgnc and len(hgnc.records) > 0:
        # probably fine to just use first match
        return hgnc.records[0].symbol
    return None


def get_gene_symbol(target_gene: TargetGene) -> str | None:
    """Acquire HGNC gene symbol given provided target gene metadata from MaveDB.

    Tokenizes the target name on whitespace and tries each token against the
    gene normalizer until one matches (gene symbols are not always the first
    token, e.g. ``"Wildtype G6PD"``).  Silently returns ``None`` if no token
    resolves — see ``_get_normalized_gene_response`` for the full description
    of this limitation.

    :param target_gene: target gene metadata given by MaveDB API
    :return: gene symbol if available
    """
    if target_gene.target_uniprot_ref:
        result = _get_hgnc_symbol(target_gene.target_uniprot_ref.id)
        if result:
            return result

    if target_gene.target_gene_name:
        for word in target_gene.target_gene_name.split(" "):
            result = _get_hgnc_symbol(word)
            if result:
                return result

    return None


def _normalize_gene(term: str) -> Gene | None:
    """Fetch normalizer response for gene term.

    :param term: gene name or referent to normalize
    :return: normalized Gene if successful
    """
    q = GeneNormalizerBuilder()
    response = q.normalize(term)
    if response.match_type > 0:
        return response.gene
    return None


def _get_normalized_gene_response(
    target_gene: TargetGene,
) -> Gene | None:
    """Fetch best normalized concept given available scoreset metadata.

    **Limitation — heuristic name parsing**: when the target name is not itself
    a valid HGNC symbol this function tokenizes it on whitespace and tries each
    token in order (e.g. ``"Wildtype G6PD"`` resolves because ``"G6PD"`` is the
    second token).  This will silently fail for names whose tokens are all
    non-HGNC strings (e.g. ``"my favourite protein"``).  When it fails,
    downstream chromosome selection has no anchor and BLAT's chromosome fallback
    may land on the wrong chromosome, causing transcript selection to return no
    results.  The only reliable fix is to ensure the target name contains the
    HGNC gene symbol as one of its whitespace-delimited tokens.

    :param target_gene: salient scoreset metadata items
    :return: Normalized gene if available
    """
    if target_gene.target_uniprot_ref:
        gene_descriptor = _normalize_gene(target_gene.target_uniprot_ref.id)
        if gene_descriptor:
            return gene_descriptor

    # Try each whitespace-delimited token from the target name. Gene symbols
    # are not always the first word (e.g. "Wildtype G6PD").
    if target_gene.target_gene_name:
        for word in target_gene.target_gene_name.split(" "):
            gene_descriptor = _normalize_gene(word)
            if gene_descriptor:
                return gene_descriptor

    return None


def _get_genomic_interval(
    extensions: list[Extension], src_name: str
) -> GeneLocation | None:
    """Extract start/end coords from extension list. Extensions in normalized genes
    can be of several different types, but we only want SequenceLocation data.

    :param extensions: extensions given in a descriptor
    :return: genomic interval if available
    """
    locations = [ext for ext in extensions if f"{src_name}_locations" in ext.name]
    if locations and len(locations[0].value) > 0:
        location_values = [
            v for v in locations[0].value if v["type"] == "SequenceLocation"
        ]
        if location_values:
            return GeneLocation(
                start=location_values[0]["start"],
                end=location_values[0]["end"],
                chromosome=get_chromosome_identifier_from_vrs_id(
                    f"ga4gh:{location_values[0]['sequenceReference']['refgetAccession']}"
                ),
            )
    return None


def get_gene_location(target_gene: TargetGene) -> GeneLocation | None:
    """Acquire gene location data from gene normalizer using metadata provided by
    scoreset.

    Delegates to ``_get_normalized_gene_response`` — see that function for the
    full description of the heuristic and its limitations.

    :param target_gene: data given by MaveDB API
    :return: gene location data if available
    """
    gene_descriptor = _get_normalized_gene_response(target_gene)
    if not gene_descriptor or not gene_descriptor.extensions:
        return None

    for src_name in ("ensembl", "ncbi"):
        loc = _get_genomic_interval(gene_descriptor.extensions, src_name)
        if loc:
            return loc

    return None


_MAX_LOCUS_INFERENCE_SPAN = 5_000_000
"""
The maximum span between the most distant variant loci of a target for which we will attempt
to infer the gene via Ensembl overlap. If the span is wider than this, we will skip locus
inference and fall back to the declared gene metadata. This threshold is set based on the
fact that even the largest human genes are well under 5 Mb in length, so a wider span likely
indicates an anomalous target (e.g. a multi-gene construct or stray positions) for which
locus inference would be unreliable and an oversized request to Ensembl would be doomed to fail.
"""


def infer_hgnc_symbol_from_genomic_loci(
    accession: str, positions: list[int]
) -> str | None:
    """Resolve the gene at a genomic accession's variant loci via Ensembl overlap.

    For a genomic-accession (``NC_``) target the variants' chromosomal positions are
    known before mapping, so the "mapper-resolved" gene can be found by asking which
    gene overlaps those loci -- the same Ensembl overlap source the post-mapping
    gene-info cascade uses, but driven from the raw positions rather than mapped VRS
    spans so it is available in time for transcript selection. The gene that selection
    stamps onto ``TxSelectResult.hgnc_symbol`` is then reused by
    ``annotate.compute_target_gene_info`` (its first-priority path), so Ensembl is
    queried once, not twice.

    Issues a **single** Ensembl overlap request over the bounding span of all loci
    (the variants of a coding target sit within one gene's span), then disambiguates
    locally using the gene coordinates the response carries -- counting how many of the
    query loci fall inside each returned gene -- so the cost is one request per target,
    not one per variant. Only ``protein_coding`` genes are considered (a coding/MANE
    transcript exists for no others). Returns the coding gene containing the most loci;
    ``None`` on a tie, no protein-coding overlap, or a bounding span too wide to be a
    single gene (see ``_MAX_LOCUS_INFERENCE_SPAN``).

    :param accession: genomic accession the variants are described against (``NC_``)
    :param positions: 0-based genomic positions of the variants
    :return: HGNC gene symbol if one gene clearly dominates, else ``None``
    """
    if not positions:
        return None

    unique_positions = sorted(set(positions))
    span = unique_positions[-1] - unique_positions[0]
    if span >= _MAX_LOCUS_INFERENCE_SPAN:
        _logger.warning(
            "Variant loci of %s span %d bp (>= %d), too wide for single-gene locus "
            "inference; falling back to declared gene metadata.",
            accession,
            span,
            _MAX_LOCUS_INFERENCE_SPAN,
        )
        return None

    try:
        features = get_overlapping_features_for_region(
            accession, unique_positions[0], unique_positions[-1] + 1, features=["gene"]
        )
    except Exception:
        _logger.exception(
            "Overlap query failed for %s over %d-%d",
            accession,
            unique_positions[0],
            unique_positions[-1],
        )
        return None

    # Disambiguate locally from the returned gene coordinates -- no extra requests.
    # Restrict to protein-coding genes: a coding transcript (MANE) only exists for one,
    # so a non-coding gene (lncRNA, pseudogene, etc.) overlapping the same locus would
    # just fail the downstream MANE lookup. Filtering here also resolves the common tie
    # between a coding gene and an overlapping non-coding one.
    counts: dict[str, int] = {}
    for feature in features:
        symbol = feature.get("external_name")
        start = feature.get("start")
        end = feature.get("end")
        biotype = feature.get("biotype")
        if not symbol or start is None or end is None or biotype != "protein_coding":
            continue

        contained = sum(1 for p in unique_positions if start <= p <= end)
        if contained:
            counts[symbol] = contained

    if not counts:
        _logger.info(
            "No protein-coding gene overlaps the variant loci of %s; cannot infer gene.",
            accession,
        )
        return None

    max_count = max(counts.values())
    candidates = sorted(g for g, c in counts.items() if c == max_count)
    if len(candidates) > 1:
        _logger.warning(
            "Multiple candidate genes overlap the variant loci of %s: %s. "
            "No gene inferred from loci.",
            accession,
            candidates,
        )
        return None

    # Ensembl's ``external_name`` is already the gene's HGNC symbol, which is what the
    # MANE-by-gene lookup matches on. Canonicalize through the gene normalizer when it
    # resolves, but fall back to the Ensembl symbol rather than dropping the gene if the
    # normalizer returns nothing -- losing a locus-confirmed gene to a normalizer miss
    # would force an avoidable "no coding transcript" skip.
    inferred = _get_hgnc_symbol(candidates[0]) or candidates[0]
    _logger.info("Inferred gene %s from variant loci of %s.", inferred, accession)
    return inferred


# --------------------------------- SeqRepo --------------------------------- #


def check_seqrepo() -> None:
    sr = get_seqrepo()
    if not sr.sr["NC_000001.11"][780000:780020]:
        msg = "SeqRepo returned no sequence for NC_000001.11 at 780000:780020. This indicates an underlying issue with SeqRepo that should be investigated."
        _logger.error(msg)
        raise DataLookupError(msg)
    conn = sr.sr.aliases._db
    try:
        # conn = sr.sr.aliases._db
        cursor = conn.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS test_table (id INTEGER PRIMARY KEY)")
        cursor.execute("INSERT INTO test_table (id) VALUES (1)")
        conn.commit()
        cursor.execute("DELETE FROM test_table WHERE id = 1")
        cursor.execute("DROP TABLE test_table")
        conn.commit()
        # conn.close()
    except sqlite3.Error as e:
        conn.close()
        _logger.error("SeqRepo sequences DB isn't writeable.")
        raise DataLookupError from e


def get_chromosome_identifier(chromosome: str) -> str:
    """Get latest NC_ accession identifier given a chromosome name.

    :param chromosome: chromosome name, e.g. ``"8"``, ``"X"``
    :return: latest ID if available
    :raise KeyError: if unable to retrieve identifier
    """
    # target sequence alignment references are chromosome names like ``"8"``, ``"X"``
    # but accession alignment information from cdot has reference accessions, beginning with "NC_"
    # for "NC_" identifiers, just return the identifier
    if chromosome.startswith("NC_"):
        return chromosome
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"
    sr = get_seqrepo()
    acs = []
    for assembly in ["GRCh38", "GRCh37"]:
        tmp_acs, _ = sr.translate_identifier(
            f"{assembly}:{chromosome}", target_namespaces="refseq"
        )
        for ac in tmp_acs:
            acs.append(ac.split("refseq:")[-1])
    if not acs:
        msg = f"Cannot retrieve NC identifier for {chromosome} from Seqrepo"
        raise KeyError(msg)

    # make sure e.g. version .10 > version .9
    sorted_results = sorted(acs, key=lambda i: int(i.split(".")[-1]))
    return sorted_results[-1]


def get_ucsc_chromosome_name(chromosome: str) -> str:
    """Get UCSC/GENCODE-style chromosome name, eg ``"chr1"`` instead of ``"1"`` or
    ``"NC_000001.11"``.

    :param chromosome: chromosome name/identifier
    :return: UCSC/GENCODE-style chromosome name
    :raise KeyError: if unable to find matching name
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(chromosome, "GRCh38")
    if not result:
        msg = (
            f"Cannot retrieve USCS-style chromosome name for {chromosome} from Seqrepo"
        )
        raise KeyError(msg)

    sorted_results = sorted([r for r in result if "chr" in r])
    try:
        return sorted_results[-1].split(":")[1]
    except IndexError as e:
        raise KeyError from e


def get_chromosome_identifier_from_vrs_id(sequence_id: str) -> str | None:
    """Get NC_ identifier given a VRS sequence ID.

    :param sequence_id: identifier a la ``ga4gh:SQ.XXXXXX``
    :return: NC_ chromosome ID
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(sequence_id, "refseq")
    if not result:
        msg = f"Cannot retrieve NC identifier for {sequence_id} from Seqrepo"
        raise KeyError(msg)

    sorted_results = sorted(result)
    return sorted_results[-1]


def get_vrs_id_from_identifier(sequence_id: str) -> str | None:
    """Get GA4GH SQ identifier given an NP_ sequence id:
    :param: GA4GH SQ digest
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(sequence_id, "ga4gh")
    if not result:
        msg = f"Cannot retrieve GA4GH SQ identifier for {sequence_id} from Seqrepo"
        raise KeyError(msg)
    sorted_results = sorted(result)
    return sorted_results[-1]


def get_sequence(
    sequence_id: str,
    start: int | None = None,
    end: int | None = None,
) -> str:
    """Get reference sequence given a sequence identifier.

    :param sequence_id: sequence identifier, e.g. ``"NP_938033.1"``
    :return: sequence
    :raise KeyError: if lookup fails
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    try:
        sequence = sr.get_sequence(sequence_id, start, end)
    except (KeyError, ValueError) as e:
        _logger.error("Unable to acquire sequence for ID: %s", sequence_id)
        raise KeyError from e
    if sequence is None:
        _logger.error("Unable to acquire sequence for ID: %s", sequence_id)
        raise KeyError
    return sequence


# -------------------------------- VRS-Python -------------------------------- #


def translate_hgvs_to_vrs(hgvs: str) -> Allele:
    """Convert HGVS variation description to VRS object.

    :param hgvs: MAVE-HGVS variation string
    :return: Corresponding VRS allele as a Pydantic class
    """
    # coerce tmp HGVS string into formally correct term
    if hgvs.startswith("NC_") and ":c." in hgvs:
        hgvs = hgvs.replace(":c.", ":g.")

    tr = TranslatorBuilder(get_seqrepo())
    allele: Allele = tr.translate_from(hgvs, "hgvs", do_normalize=False)

    if (
        not isinstance(allele.location, SequenceLocation)
        or not isinstance(allele.location.start, int)
        or not isinstance(allele.location.end, int)
        or not isinstance(allele.state, LiteralSequenceExpression)
    ):
        raise ValueError
    return allele


def build_ref_identical_allele(sequence_id: str) -> Allele:
    """Build a whole-sequence reference-identical VRS Allele from a sequence identifier.

    Accepts either a GA4GH SQ digest (``SQ.xxx``, without the ``ga4gh:`` prefix) or a
    named accession such as ``NP_``, ``NM_``, or ``NC_``.

    :param sequence_id: GA4GH SQ digest or named accession
    :return: VRS Allele spanning the full sequence with a ReferenceLengthExpression state
    :raises DataLookupError: if the sequence identifier or metadata lookup fails
    """
    if sequence_id.startswith("SQ."):
        ga4gh_id = f"ga4gh:{sequence_id}"
    else:
        try:
            ga4gh_id = get_vrs_id_from_identifier(sequence_id)
        except KeyError as e:
            msg = f"Could not retrieve GA4GH identifier for accession {sequence_id!r}"
            _logger.error(msg)
            raise DataLookupError(msg) from e

    sr = get_seqrepo()
    try:
        metadata = sr.get_metadata(ga4gh_id)
    except KeyError as e:
        msg = f"Could not retrieve metadata for sequence {ga4gh_id!r}"
        _logger.error(msg)
        raise DataLookupError(msg) from e

    length = metadata["length"]
    refget_accession = ga4gh_id.split("ga4gh:")[-1]

    seq_ref = SequenceReference(refgetAccession=refget_accession)
    location = SequenceLocation(sequenceReference=seq_ref, start=0, end=length)
    state = ReferenceLengthExpression(length=length, repeatSubunitLength=length)

    return Allele(location=location, state=state)


def translate_ref_identical_to_vrs(hgvs_string: str) -> Allele:
    """Convert a reference-identical HGVS variant to a VRS Allele.

    Handles reference-identical variants such as ``NM_001234.1:c.=``,
    ``NP_001234.1:p.=``, and ``NC_000001.11:g.=``, which regular VRS
    translation does not support. Returns an Allele with a
    ``ReferenceLengthExpression`` state spanning the full reference sequence.

    :param hgvs_string: HGVS reference-identical variant string (e.g. ``NM_001234.1:c.=``)
    :return: VRS Allele spanning the full reference sequence
    :raises ValueError: if ``hgvs_string`` is not a valid reference-identical HGVS expression
    :raises DataLookupError: if the sequence identifier or metadata lookup fails
    """
    if ":" not in hgvs_string or not hgvs_string.endswith(".="):
        msg = f"Not a reference-identical HGVS expression: {hgvs_string!r}"
        raise ValueError(msg)

    accession = hgvs_string.split(":")[0]
    return build_ref_identical_allele(accession)


# ------------------------------ HGVS Projection ------------------------------ #


def project_genomic_hgvs_to_coding(g_hgvs: str, transcript: str) -> str:
    """Project a genomic HGVS expression onto a coding transcript (``g.`` -> ``c.``).

    Uses the cdot-backed hgvs ``VariantMapper`` already wired up for this package
    (see :func:`init_hgvs_tools`). This is the mapper's *projection* of a measured
    variant onto its own coding form -- distinct from the equivalence-class
    expansion (all synonymous codons) that the reverse-translation job owns.

    :param g_hgvs: genomic HGVS string, e.g. ``NC_000001.11:g.123A>G``
    :param transcript: coding transcript accession to project onto, e.g. ``NM_...``
    :return: coding HGVS string (``NM_...:c....``)
    """
    tools = TranslatorBuilder(get_seqrepo()).hgvs_tools
    var_g = tools.parser.parse(g_hgvs)
    var_c = tools.variant_mapper.g_to_c(var_g, transcript)
    return str(var_c)


def coding_hgvs_is_intronic(c_hgvs: str) -> bool:
    """Return whether a coding HGVS expression refers to an intronic position.

    Detected from the parsed base-offset position (a non-zero intron offset), not the
    string: a textual ``-`` would misclassify 5'UTR positions (``c.-20``, offset 0) as
    intronic. Intronic variants have no protein consequence and cannot be represented by
    ga4gh's VRS tools, so callers projecting a variant onto a transcript skip them.

    :param c_hgvs: coding HGVS string, e.g. ``NM_...:c.2002-1_2003del``
    :return: ``True`` if either endpoint carries a non-zero intron offset
    """
    tools = TranslatorBuilder(get_seqrepo()).hgvs_tools
    pos = tools.parser.parse(c_hgvs).posedit.pos
    return bool(getattr(pos.start, "offset", 0) or getattr(pos.end, "offset", 0))


def project_coding_hgvs_to_protein(c_hgvs: str) -> str:
    """Project a coding HGVS expression onto its protein consequence (``c.`` -> ``p.``).

    :param c_hgvs: coding HGVS string, e.g. ``NM_...:c.76A>G``
    :return: protein HGVS string (``NP_...:p....``)
    """
    tools = TranslatorBuilder(get_seqrepo()).hgvs_tools
    var_c = tools.parser.parse(c_hgvs)
    var_p = tools.variant_mapper.c_to_p(var_c)
    return str(var_p)


def get_genomic_accession_for_transcript(
    transcript: str, assembly: str = "GRCh38"
) -> str | None:
    """Resolve the genomic contig (``NC_``) a coding transcript aligns to on an assembly.

    Needed for the ``c.`` -> ``g.`` projection: unlike :func:`project_genomic_hgvs_to_coding`
    (whose ``g_to_c`` reads the genomic accession from the input string), hgvs's
    ``c_to_g`` requires the target genomic accession passed explicitly. cdot carries one
    contig per genome build in the transcript record, but its public
    ``get_tx_mapping_options`` flattens the build keys away, so disambiguate by
    intersecting the candidate contigs with the requested assembly's contig set
    (``get_assembly_map``).

    :param transcript: coding transcript accession, e.g. ``NM_004333.6``
    :param assembly: genome build to resolve the contig on (``GRCh38`` or ``GRCh37``)
    :return: the assembly's ``NC_`` contig for the transcript, or ``None`` if unresolvable
    """
    cd = cdot_rest()
    try:
        options = cd.get_tx_mapping_options(transcript)
        assembly_contigs = cd.get_assembly_map(assembly)
    except Exception:
        _logger.exception(
            "Could not resolve genomic contig for transcript %s on %s",
            transcript,
            assembly,
        )
        return None

    for option in options:
        if option["alt_ac"] in assembly_contigs:
            return option["alt_ac"]

    return None


def project_coding_hgvs_to_genomic(c_hgvs: str, alt_ac: str) -> str:
    """Project a coding HGVS expression onto its genomic form (``c.`` -> ``g.``).

    The deterministic genomic re-expression of a coding variant, on the contig
    ``alt_ac`` (resolve it with :func:`get_genomic_accession_for_transcript`). Like the
    other projection helpers this routes through the cdot-backed package hgvs
    ``VariantMapper`` (see :func:`init_hgvs_tools`); ``c_to_g`` needs the target contig
    explicitly because the coding expression does not name one.

    :param c_hgvs: coding HGVS string, e.g. ``NM_004333.6:c.1A>G``
    :param alt_ac: genomic contig accession to project onto, e.g. ``NC_000007.14``
    :return: genomic HGVS string (``NC_...:g....``)
    """
    tools = TranslatorBuilder(get_seqrepo()).hgvs_tools
    var_c = tools.parser.parse(c_hgvs)
    var_g = tools.variant_mapper.c_to_g(var_c, alt_ac)
    return str(var_g)


# ----------------------------------- MANE ----------------------------------- #


def _sort_mane_result(description: ManeDescription) -> int:
    if description.transcript_priority == TranscriptPriority.MANE_SELECT:
        return 2
    if description.transcript_priority == TranscriptPriority.MANE_PLUS_CLINICAL:
        return 1

    # should be unreachable.
    _logger.warning(
        "Unrecognized transcript priority value %s for transcript description of %s",
        description.transcript_priority,
        description.refseq_nuc,
    )
    return 0


def _mane_row_to_description(row: dict[str, Any]) -> ManeDescription:
    """Build a ``ManeDescription`` from a row of the MANE summary dataframe.

    Both the transcript-keyed and gene-keyed lookups read the same dataframe
    (``ManeTranscriptMappings.df``), so they share column names.
    """
    return ManeDescription(
        ncbi_gene_id=row["#NCBI_GeneID"],
        ensembl_gene_id=row["Ensembl_Gene"],
        hgnc_gene_id=row["HGNC_ID"],
        symbol=row["symbol"],
        name=row["name"],
        refseq_nuc=row["RefSeq_nuc"],
        refseq_prot=row["RefSeq_prot"],
        ensembl_nuc=row["Ensembl_nuc"],
        ensembl_prot=row["Ensembl_prot"],
        transcript_priority=TranscriptPriority(
            "_".join(row["MANE_status"].lower().split())
        ),
        grch38_chr=row["GRCh38_chr"],
        chr_start=row["chr_start"],
        chr_end=row["chr_end"],
        chr_strand=row["chr_strand"],
    )


def get_mane_transcripts(transcripts: list[str]) -> list[ManeDescription]:
    """Get corresponding MANE data for transcripts. Results given in order of
    transcript preference.

    :param transcripts: candidate transcripts list
    :return: complete MANE descriptions
    """
    mane_df = CoolSeqToolBuilder().mane_transcript_mappings.df
    mane_results = mane_df.filter(pl.col("RefSeq_nuc").is_in(transcripts))
    mane_data = [_mane_row_to_description(row) for row in mane_results.rows(named=True)]
    mane_data.sort(key=_sort_mane_result)
    return mane_data


def get_mane_transcripts_for_gene(gene_symbol: str) -> list[ManeDescription]:
    """Get MANE transcript(s) for a gene symbol directly, without a prior alignment.

    Uses cool-seq-tool's gene-keyed MANE lookup. Unlike ``get_mane_transcripts``
    (which filters by candidate transcript accessions seeded from a genomic
    alignment), this resolves the gene's MANE transcript straight from the gene
    symbol -- the path needed for accession-based genomic (``NC_``) targets that
    have no BLAT alignment to seed candidate transcripts.

    :param gene_symbol: HGNC gene symbol
    :return: MANE descriptions, MANE Select first then MANE Plus Clinical
    """
    rows = CoolSeqToolBuilder().mane_transcript_mappings.get_gene_mane_data(gene_symbol)
    mane_data = [_mane_row_to_description(row) for row in rows]
    mane_data.sort(key=_sort_mane_result)
    return mane_data


# --------------------------------- Ensembl --------------------------------- #


def get_overlapping_features_for_region(
    chromosome: str, start: int, end: int, features: list[str] | None = None
) -> list[dict[str, Any]]:
    """Get genes overlapping a specific genomic region.

    :param chromosome: Chromosome identifier
    :param start: Start position of the region
    :param end: End position of the region
    :param features: List of features to retrieve (default is ["gene"])
    :return: List of overlapping gene symbols
    """
    if not features:
        features = ["gene"]
        _logger.debug("No features specified, defaulting to %s", features)

    chrom = get_chromosome_identifier(chromosome)

    query = f"/{chrom}:{start}-{end}"
    if features:
        query += "?"
    for feature in features:
        query += f"feature={feature};"

    try:
        _logger.debug(
            "Fetching overlapping features for region %s:%d-%d with features %s",
            chromosome,
            start,
            end,
            features,
        )

        url = f"{ENSEMBL_API_URL}/overlap/region/human{query}"
        response = request_with_backoff(
            url, headers={"Content-Type": "application/json"}
        )
        response.raise_for_status()
    except httpx.HTTPError as e:
        _logger.error(
            "Failed to fetch overlapping features for region %s-%s on chromosome %s: %s",
            start,
            end,
            chromosome,
            e,
        )
        return []

    overlapping_features = response.json()
    _logger.debug(
        "Successfully fetched %d overlapping features for region %s:%d-%d with features %s",
        len(overlapping_features),
        chromosome,
        start,
        end,
        features,
    )
    return overlapping_features


# ---------------------------------- Misc. ---------------------------------- #


def get_uniprot_sequence(uniprot_id: str) -> str | None:
    """Get sequence directly from UniProt.

    :param uniprot_id: ID provided with target info
    :return: transcript accession if successful
    :raise HTTPError: if response comes with an HTTP error code
    """
    url = f"https://www.ebi.ac.uk/proteins/api/proteins?accession={uniprot_id.split(':')[1]}&format=json"
    response = httpx.get(url, timeout=30)
    response.raise_for_status()
    json = response.json()
    return json[0]["sequence"]["sequence"]
