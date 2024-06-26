"""Map transcripts to VRS objects."""

import logging
from collections.abc import Iterable
from itertools import cycle

import click
from Bio.Seq import Seq
from bioutils.accessions import infer_namespace
from cool_seq_tool.schemas import AnnotationLayer, Strand
from ga4gh.core import ga4gh_identify, sha512t24u
from ga4gh.vrs._internal.models import (
    Allele,
    Haplotype,
    LiteralSequenceExpression,
    ReferenceLengthExpression,
    SequenceLocation,
    SequenceString,
)
from ga4gh.vrs.normalize import normalize
from mavehgvs.util import parse_variant_strings
from mavehgvs.variant import Variant

from dcd_mapping.lookup import (
    cdot_rest,
    get_chromosome_identifier,
    get_seqrepo,
    translate_hgvs_to_vrs,
)
from dcd_mapping.schemas import (
    AlignmentResult,
    MappedScore,
    ScoreRow,
    ScoresetMetadata,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
)

__all__ = ["vrs_map", "VrsMapError"]


_logger = logging.getLogger(__name__)


class VrsMapError(Exception):
    """Raise in case of VRS mapping errors."""


def _hgvs_variant_is_valid(hgvs_string: str) -> bool:
    return not hgvs_string.endswith((".=", ")", "X"))


def _process_any_aa_code(hgvs_pro_string: str) -> str:
    """Substitute "Xaa" for "?" in variation expression.

    Some expressions seem to use the single-character "?" wildcard in the context of
    three-letter amino acid codes. This is weird, and the proper replacement is "Xaa".

    Note that we currently do NOT make any alterations to nucleotide strings that use
    weird apparently-wildcard characters like "X" -- we just treat them as invalid (see
    _hgvs_variant_is_valid()).

    :param hgvs_string: MAVE HGVS expression
    :return: processed variation (equivalent to input if no wildcard code found)
    """
    if "?" in hgvs_pro_string:
        _logger.debug("Substituting Xaa for ? in %s", hgvs_pro_string)
        hgvs_pro_string = hgvs_pro_string.replace("?", "Xaa")
    return hgvs_pro_string


def _create_pre_mapped_hgvs_strings(
    raw_description: str,
    layer: AnnotationLayer,
    tx: TxSelectResult | None = None,
    alignment: AlignmentResult | None = None,
    accession_id: str | None = None,
) -> list[str]:
    """Generate a list of (pre-mapped) HGVS strings from one long string containing many valid HGVS substrings

    Currently, the provided transcript is used as the reference for the hgvs string, but this is inaccurate
    because pre-mapped variants should be relative to the user-provided target sequence, not an external accession.
    Any offset between the transcript and target sequence is not taken into account here (the variant position
    is relative to the target sequence).

    :param raw_description: A string containing valid HGVS sub-strings
    :param layer: An enum denoting the targeted annotation layer of these HGVS strings
    :param tx: A TxSelectResult object defining the transcript we are mapping to (or None).
    :param alignment: An AlignmentResult object defining the alignment we are mapping to (or None).
    :return: A list of HGVS strings prior to being mapped to the `tx` or `alignment`
    """
    if layer is AnnotationLayer.PROTEIN and tx is None:
        msg = f"Transcript result must be provided for {layer} annotations."
        raise ValueError(msg)
    if layer is AnnotationLayer.GENOMIC and alignment is None and accession_id is None:
        msg = f"Alignment result or accession id must be provided for {layer} annotations."
        raise ValueError(msg)

    raw_variant_strings = _parse_raw_variant_str(raw_description)
    variants, errors = parse_variant_strings(raw_variant_strings)

    hgvs_strings = []
    for variant, error in zip(variants, errors, strict=True):
        if error is not None:
            msg = f"Variant could not be parsed by mavehgvs: {error}"
            raise ValueError(msg)

        if accession_id:
            hgvs_strings.append(accession_id + ":" + str(variant))
        # Ideally we would create an HGVS string namespaced to GA4GH. The line below
        # creates such a string, but it is not able to be parsed by the GA4GH VRS translator.
        # hgvs_strings.append('ga4gh:' + sequence_id + ':' + str(variant))
        elif layer is AnnotationLayer.PROTEIN:
            assert tx  # noqa: S101. mypy help
            hgvs_strings.append(tx.np + ":" + str(variant))
        elif layer is AnnotationLayer.GENOMIC:
            assert alignment  # noqa: S101. mypy help
            hgvs_strings.append(
                get_chromosome_identifier(alignment.chrom) + ":" + str(variant)
            )
        else:
            msg = (
                f"Could not generate HGVS strings for invalid AnnotationLayer: {layer}"
            )
            raise ValueError(msg)

    return hgvs_strings


def _create_post_mapped_hgvs_strings(
    raw_description: str,
    layer: AnnotationLayer,
    tx: TxSelectResult | None = None,
    alignment: AlignmentResult | None = None,
) -> list[str]:
    """Generate a list of (post-mapped) HGVS strings from one long string containing many valid HGVS substrings.

    For protein annotations, these strings must be adjusted to match the offset defined by the start of the
    transcript sequence. For genomic annotations, these strings must be adjusted to match the coordinates of
    the reference alignment.

    :param raw_description: A string containing valid HGVS sub-strings
    :param layer: An enum denoting the targeted annotation layer of these HGVS strings
    :param tx: A TxSelectResult object defining the transcript we are mapping to (or None)
    :param alignment: An AlignmentResult object defining the alignment we are mapping to (or None)
    :return: A list of HGVS strings relative to the `tx` or `alignment`
    """
    if layer is AnnotationLayer.PROTEIN and tx is None:
        msg = f"Transcript result must be provided for {layer} annotations (Transcript was `{tx}`)."
        raise ValueError(msg)
    if layer is AnnotationLayer.GENOMIC and alignment is None:
        msg = f"Alignment result must be provided for {layer} annotations (Alignment was `{alignment}`)."
        raise ValueError(msg)

    raw_variants = _parse_raw_variant_str(raw_description)
    variants, errors = parse_variant_strings(raw_variants)

    hgvs_strings = []
    for variant, error in zip(variants, errors, strict=True):
        if error is not None:
            msg = f"Variant could not be parsed by mavehgvs: {error}"
            raise ValueError(msg)

        if layer is AnnotationLayer.PROTEIN:
            assert tx  # noqa: S101. mypy help

            variant = _adjust_protein_variant_to_ref(variant, tx)
            hgvs_strings.append(tx.np + ":" + str(variant))
        elif layer is AnnotationLayer.GENOMIC:
            assert alignment  # noqa: S101. mypy help

            variant = _adjust_genomic_variant_to_ref(variant, alignment)
            hgvs_strings.append(
                get_chromosome_identifier(alignment.chrom) + ":" + str(variant)
            )
        else:
            msg = (
                f"Could not generate HGVS strings for invalid AnnotationLayer: {layer}"
            )
            raise ValueError(msg)

    return hgvs_strings


def _adjust_protein_variant_to_ref(
    variant: Variant,
    tx: TxSelectResult,
) -> Variant:
    if isinstance(variant.positions, Iterable):
        for position in variant.positions:
            position.position = position.position + tx.start
        return variant

    variant.positions.position = variant.positions.position + tx.start
    return variant


def _adjust_genomic_variant_to_ref(
    variant: Variant,
    alignment: AlignmentResult,
) -> Variant:
    """Adjust a variant relative to a provided alignment.

    :param variant: A variant object relative to a scoreset's target sequence
    :param alignment: An AlignmentResult object denoting the alignment we are mapping to
    :return: A variant object that describes the variant relative to the provided alignment result
    """
    # adjust starts - hgvs uses 1-based numbering for c. sequences, while blat hits are 0-based
    starts = []
    if isinstance(variant.positions, Iterable):
        is_multi_position = True
        for position in variant.positions:
            starts.append(position.position - 1)
    else:
        is_multi_position = False
        starts.append(variant.positions.position - 1)

    # get hit
    query_subrange_containing_hit = None
    target_subrange_containing_hit = None
    for query_subrange, target_subrange in zip(
        alignment.query_subranges, alignment.hit_subranges, strict=True
    ):
        if all(
            start >= query_subrange.start and start < query_subrange.end
            for start in starts
        ):
            query_subrange_containing_hit = query_subrange
            target_subrange_containing_hit = target_subrange
            break

    if query_subrange_containing_hit is None or target_subrange_containing_hit is None:
        msg = "Variant was not contained, or multi-position variant was not fully contained, within the aligned portion of the query sequence."
        raise ValueError(msg)

    for idx, start in enumerate(starts):
        if alignment.strand is Strand.POSITIVE:
            # get variant start relative to the reference (the "hit")
            # distance from beginning of query to variant start position:
            query_to_start = start - query_subrange_containing_hit.start

            # distance from beginning of ref to the variant start position:
            ref_to_start = target_subrange_containing_hit.start + query_to_start
        else:
            # picture the rev comp of the query/variant as mapping to the positive strand of the ref
            # the start of the reverse complement of the variant is the end of the "original" variant
            # so we need to know where the end of the original variant is, relative to the query molecule
            end = start

            # subtract 1 from end of hit range, because blat ranges are 0-based [start, end)
            ref_to_start = (target_subrange_containing_hit.end - 1) - (
                end - query_subrange_containing_hit.start
            )

        # add distance from ref to variant start; hgvs is 1-based, so convert back to 1-based
        if is_multi_position:
            variant.positions[idx].position = ref_to_start + 1
        else:
            variant.positions.position = ref_to_start + 1

    # get reverse complement of sequence if the target maps to the negative strand of the reference
    if alignment.strand is Strand.NEGATIVE:
        # variant._sequences can be a string or an iterable
        if isinstance(variant._sequences, str):
            variant._sequences = str(Seq(variant._sequences).reverse_complement())
        elif variant._sequences is not None:
            revcomp_sequences_list = []
            for sequence in variant._sequences:
                revcomp_sequences_list.append(str(Seq(sequence).reverse_complement()))
            variant._sequences = revcomp_sequences_list

        # reverse order of positions tuple
        if is_multi_position:
            variant._positions = tuple(reversed(list(variant.positions)))

    # change prefix from c. to g. since variant is now relative to chr reference
    variant._prefix = "g"

    return variant


def _parse_raw_variant_str(raw_description: str) -> list[str]:
    """Parse a string which may contain many HGVS strings into a list of each one.

    :param raw_description: A string that may contain a list of variant descriptions or a single variant description
    :return: A list of HGVS strings
    """
    if "[" in raw_description:
        prefix = raw_description[0:2]
        return [prefix + var for var in set(raw_description[3:-1].split(";"))]

    return [raw_description]


def _map_protein_coding_pro(
    row: ScoreRow,
    sequence_id: str,
    transcript: TxSelectResult,
) -> MappedScore | None:
    """Construct VRS object mapping for ``hgvs_pro`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row: A row of output from a MaveDB score set
    :param sequence_id: The GA4GH accession for the provided sequence
    :param transcript: The transcript selection information for a score set
    :return: VRS mapping object if mapping succeeds
    """
    if (
        row.hgvs_pro in {"_wt", "_sy", "NA"}
        or "fs" in row.hgvs_pro
        or len(row.hgvs_pro) == 3
    ):
        _logger.warning(
            "Can't process variant syntax %s for %s", row.hgvs_pro, row.accession
        )
        return None

    # TODO: Handle edge cases without hardcoding URNs.
    # Special case for experiment set urn:mavedb:0000097
    if row.hgvs_pro.startswith("NP_009225.1:p."):
        vrs_variation = translate_hgvs_to_vrs(row.hgvs_pro)
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            annotation_layer=AnnotationLayer.PROTEIN,
            pre_mapped=vrs_variation,
            post_mapped=vrs_variation,
        )

    pre_mapped_hgvs_strings = _create_pre_mapped_hgvs_strings(
        row.hgvs_pro,
        AnnotationLayer.PROTEIN,
        tx=transcript,
    )
    post_mapped_hgvs_strings = _create_post_mapped_hgvs_strings(
        row.hgvs_pro,
        AnnotationLayer.PROTEIN,
        tx=transcript,
    )

    pre_mapped_protein = _construct_vrs_allele(
        pre_mapped_hgvs_strings,
        AnnotationLayer.PROTEIN,
        sequence_id,
        True,
    )
    post_mapped_protein = _construct_vrs_allele(
        post_mapped_hgvs_strings,
        AnnotationLayer.PROTEIN,
        None,
        False,
    )

    if pre_mapped_protein and post_mapped_protein:
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            annotation_layer=AnnotationLayer.PROTEIN,
            pre_mapped=pre_mapped_protein,
            post_mapped=post_mapped_protein,
        )

    return None


def _map_genomic(
    row: ScoreRow,
    sequence_id: str,
    align_result: AlignmentResult | None,
) -> MappedScore | None:
    """Construct VRS object mapping for ``hgvs_nt`` variant column entry

    These arguments are a little lazy and could be pruned down later

    :param row: A row of output from a MaveDB score set
    :param sequence_id: The GA4GH accession (if target-based score set), or RefSeq/Ensembl accession associated with score set
    :param align_result: The transcript selection information for a score set
    :return: VRS mapping object if mapping succeeds
    """
    namespace = infer_namespace(sequence_id).lower()
    if align_result is None:
        # for contig accession based score sets, no mapping is performed,
        # so pre- and post-mapped alleles are the same
        pre_mapped_hgvs_strings = (
            post_mapped_hgvs_strings
        ) = _create_pre_mapped_hgvs_strings(
            row.hgvs_nt,
            AnnotationLayer.GENOMIC,
            accession_id=sequence_id,
        )

    elif namespace in ("refseq", "ncbi", "ensembl"):
        # nm/enst way
        pre_mapped_hgvs_strings = _create_pre_mapped_hgvs_strings(
            row.hgvs_nt,
            AnnotationLayer.GENOMIC,
            accession_id=sequence_id,
        )
        post_mapped_hgvs_strings = _create_post_mapped_hgvs_strings(
            row.hgvs_nt,
            AnnotationLayer.GENOMIC,
            alignment=align_result,
        )
    elif namespace == "ga4gh":
        # target seq way
        pre_mapped_hgvs_strings = _create_pre_mapped_hgvs_strings(
            row.hgvs_nt,
            AnnotationLayer.GENOMIC,
            alignment=align_result,
        )
        post_mapped_hgvs_strings = _create_post_mapped_hgvs_strings(
            row.hgvs_nt,
            AnnotationLayer.GENOMIC,
            alignment=align_result,
        )
    else:
        msg = f"Namespace not supported: {namespace}"
        raise ValueError(msg)

    pre_mapped_genomic = _construct_vrs_allele(
        pre_mapped_hgvs_strings,
        AnnotationLayer.GENOMIC,
        sequence_id,
        True,
    )
    post_mapped_genomic = _construct_vrs_allele(
        post_mapped_hgvs_strings,
        AnnotationLayer.GENOMIC,
        None,
        False,
    )

    if pre_mapped_genomic and post_mapped_genomic:
        return MappedScore(
            accession_id=row.accession,
            score=row.score,
            annotation_layer=AnnotationLayer.GENOMIC,
            pre_mapped=pre_mapped_genomic,
            post_mapped=post_mapped_genomic,
        )

    return None


def _get_allele_sequence(allele: Allele) -> str:
    """Get sequence for Allele

    :param allele: VRS allele
    :return: sequence
    :raise ValueError: if sequence is none
    """
    dp = get_seqrepo()
    start = allele.location.start
    end = allele.location.end
    sequence = dp.get_sequence(
        f"ga4gh:{allele.location.sequenceReference.refgetAccession}", start, end
    )
    if sequence is None:
        raise ValueError
    return sequence


def store_sequence(sequence: str) -> str:
    """Store sequence in SeqRepo.

    :param sequence: raw sequence (ie nucleotides or amino acids)
    :return: sequence ID (sans prefix, which is ``"ga4gh"``)
    """
    sequence_id = f"SQ.{sha512t24u(sequence.encode('ascii'))}"
    alias_dict_list = [{"namespace": "ga4gh", "alias": sequence_id}]
    sr = get_seqrepo()
    sr.sr.store(sequence, alias_dict_list)
    return sequence_id


def _hgvs_nt_is_valid(hgvs_nt: str) -> bool:
    """Check for invalid or unavailable nucleotide MAVE-HGVS variation

    :param hgvs_nt: MAVE_HGVS nucleotide expression
    :return: True if expression appears populated and valid
    """
    return (
        (hgvs_nt != "NA")
        and (hgvs_nt not in {"_wt", "_sy", "="})
        and (len(hgvs_nt) != 3)
    )


def _map_protein_coding(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    transcript: TxSelectResult,
    align_result: AlignmentResult,
) -> list[MappedScore]:
    """Perform mapping on protein coding experiment results

    :param metadata: The metadata for a score set
    :param records: The list of MAVE variants in a given score set
    :param transcript: The transcript data for a score set
    :param align_results: The alignment data for a score set
    :return: A list of mappings
    """
    if metadata.target_sequence_type == TargetSequenceType.DNA:
        sequence = str(
            Seq(metadata.target_sequence).translate(table="1", stop_symbol="")
        )
        psequence_id = store_sequence(sequence)
        gsequence_id = store_sequence(metadata.target_sequence)
    else:
        sequence = metadata.target_sequence
        psequence_id = gsequence_id = store_sequence(sequence)

    variations: list[MappedScore] = []
    for row in records:
        hgvs_pro_mappings = _map_protein_coding_pro(row, psequence_id, transcript)
        if hgvs_pro_mappings:
            variations.append(hgvs_pro_mappings)
        else:
            _logger.warning(
                "Encountered apparently invalid protein variants in %s: %s",
                row.accession,
                row.hgvs_pro,
            )

        if _hgvs_nt_is_valid(row.hgvs_nt):
            hgvs_nt_mappings = _map_genomic(row, gsequence_id, align_result)

            if hgvs_nt_mappings:
                variations.append(hgvs_nt_mappings)
            else:
                _logger.warning(
                    "Encountered apparently invalid genomic variants in %s: %s",
                    row.accession,
                    row.hgvs_nt,
                )

    return variations


def _map_regulatory_noncoding(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    align_result: AlignmentResult,
) -> list[MappedScore]:
    """Perform mapping on noncoding/regulatory experiment results

    :param metadata: metadata for URN
    :param records: list of MAVE experiment result rows
    :param align_result: An AlignmentResult object for a score set
    :return: A list of VRS mappings
    """
    variations: list[MappedScore] = []
    sequence_id = store_sequence(metadata.target_sequence)

    for row in records:
        if (
            row.hgvs_nt in {"_wt", "_sy", "="}
            or "fs" in row.hgvs_nt
            or len(row.hgvs_nt) == 3
        ):
            _logger.warning(
                "Can't process variant syntax %s for %s", row.hgvs_nt, metadata.urn
            )
            continue

        hgvs_nt_mappings = _map_genomic(row, sequence_id, align_result)

        if hgvs_nt_mappings:
            variations.append(hgvs_nt_mappings)
        else:
            _logger.warning(
                "Encountered apparently invalid genomic variants in %s: %s",
                row.accession,
                row.hgvs_nt,
            )

    return variations


def store_accession(
    accession_id: str,
) -> None:
    namespace = infer_namespace(accession_id)
    alias_dict_list = [{"namespace": namespace, "alias": accession_id}]
    cd = cdot_rest()
    sequence = cd.get_seq(accession_id)
    sr = get_seqrepo()
    sr.sr.store(sequence, alias_dict_list)


def _map_accession(
    metadata: ScoresetMetadata,
    records: list[ScoreRow],
    align_result: AlignmentResult,
) -> list[MappedScore]:
    variations: list[MappedScore] = []
    # see if accession is in seqrepo
    # if not, get seq from cdot and add to seqrepo
    sequence_id = metadata.target_accession
    if sequence_id is None:
        raise ValueError

    store_accession(sequence_id)

    for row in records:
        if (
            row.hgvs_nt in {"_wt", "_sy", "="}
            or "fs" in row.hgvs_nt
            or len(row.hgvs_nt) == 3
        ):
            _logger.warning(
                "Can't process variant syntax %s for %s", row.hgvs_nt, metadata.urn
            )
            continue

        hgvs_nt_mappings = _map_genomic(row, sequence_id, align_result)

        if hgvs_nt_mappings:
            variations.append(hgvs_nt_mappings)
        else:
            _logger.warning(
                "Encountered apparently invalid genomic variants in %s: %s",
                row.accession,
                row.hgvs_nt,
            )

    return variations


def _rle_to_lse(
    rle: ReferenceLengthExpression, location: SequenceLocation
) -> LiteralSequenceExpression:
    """Coerce ReferenceLengthExpression to LiteralSequenceExpression.

    RLEs are helpful for long repeating sequences but a) unnecessary here and b)
    create incompatibilities with some data extraction further down so to simplify,
    we'll just turn them into equivalent LiteralSequenceExpressions.
    """
    sr = get_seqrepo()
    sequence_id = location.sequenceReference.refgetAccession
    start: int = location.start
    end = start + rle.repeatSubunitLength
    subsequence = sr.get_sequence(f"ga4gh:{sequence_id}", start, end)
    c = cycle(subsequence)
    derived_sequence = ""
    for _ in range(rle.length):
        derived_sequence += next(c)
    return LiteralSequenceExpression(sequence=derived_sequence)


def _construct_vrs_allele(
    hgvs_strings: list[str],
    layer: AnnotationLayer,
    sequence_id: str | None,
    pre_map: bool,
) -> Allele | Haplotype:
    alleles: list[Allele] = []
    for hgvs_string in hgvs_strings:
        allele = translate_hgvs_to_vrs(hgvs_string)

        if pre_map:
            if sequence_id is None:
                msg = "Must provide sequence id to construct pre-mapped VRS allele"
                raise ValueError(msg)
            allele.location.sequenceReference.refgetAccession = sequence_id

        if "dup" in hgvs_string:
            allele.state.sequence = SequenceString(2 * _get_allele_sequence(allele))

        # TODO check assumption that c.= leads to an "N" in the sequence.root
        if allele.state.sequence.root == "N" and layer == AnnotationLayer.GENOMIC:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele))

        if "=" in hgvs_string and layer == AnnotationLayer.PROTEIN:
            allele.state.sequence = SequenceString(_get_allele_sequence(allele))

        allele = normalize(allele, data_proxy=get_seqrepo())

        if isinstance(allele.state, ReferenceLengthExpression):
            _logger.debug(
                "Coercing state for %s into LSE: %s",
                hgvs_string,
                allele.state.model_dump_json(),
            )
            allele.state = _rle_to_lse(allele.state, allele.location)

        # Run ga4gh_identify to assign VA digest
        allele.id = ga4gh_identify(allele)
        alleles.append(allele)

    if not alleles:
        msg = f"Input variant hgvs_string(s) could not be translated to an allele: {hgvs_strings}."
        raise ValueError(msg)

    if len(alleles) > 1:
        return Haplotype(members=alleles)

    return alleles[0]


def vrs_map(
    metadata: ScoresetMetadata,
    align_result: AlignmentResult,
    records: list[ScoreRow],
    transcript: TxSelectResult | None = None,
    silent: bool = True,
) -> list[MappedScore] | None:
    """Given a description of a MAVE scoreset and an aligned transcript, generate
    the corresponding VRS objects.

    :param metadata: salient MAVE scoreset metadata
    :param align_result: output from the sequence alignment process
    :param records: scoreset records
    :param transcript: output of transcript selection process
    :param silent: If true, suppress console output
    :return: A list of mapping results
    """
    if metadata.urn == "urn:mavedb:00000072-a-1":
        msg = f"No RefSeq accession is available for {metadata.urn}."
        if not silent:
            click.echo(msg)
        _logger.warning(msg)
        return None

    if metadata.target_accession:
        return _map_accession(metadata, records, align_result)

    if metadata.target_gene_category == TargetType.PROTEIN_CODING and transcript:
        return _map_protein_coding(
            metadata,
            records,
            transcript,
            align_result,
        )
    return _map_regulatory_noncoding(
        metadata,
        records,
        align_result,
    )
