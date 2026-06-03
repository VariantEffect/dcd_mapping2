"""Tests for VRS allele identification helpers.

These tests run offline: ``ga4gh_identify`` is pure computation given a
well-formed VRS object, and ``normalize`` is patched so no SeqRepo or UTA
access is required.
"""

import pytest
from ga4gh.vrs._internal.models import (
    Allele,
    LiteralSequenceExpression,
    SequenceLocation,
    SequenceReference,
)

from dcd_mapping.vrs_utils import identify_allele, normalize_and_identify


def _make_allele(start: int = 1, end: int = 2, sequence: str = "A") -> Allele:
    """Build a minimal in-memory VRS Allele with no cached digests."""
    return Allele(
        location=SequenceLocation(
            sequenceReference=SequenceReference(
                refgetAccession="SQ.0123456789abcdef0123456789abcdef"
            ),
            start=start,
            end=end,
        ),
        state=LiteralSequenceExpression(sequence=sequence),
    )


def test_identify_allele_returns_va_digest():
    assert identify_allele(_make_allele()).startswith("ga4gh:VA.")


def test_identify_allele_is_content_addressed():
    # Identical content -> identical digest; differing content -> differing digest.
    a = _make_allele()
    b = _make_allele()
    c = _make_allele(start=5, end=6)
    assert identify_allele(a) == identify_allele(b)
    assert identify_allele(a) != identify_allele(c)


def test_identify_allele_clears_stale_digests():
    """The whole point of the helper: a cached pre-mutation digest must not
    leak into the post-mutation identifier.

    Reproduces the bug pattern where an allele is constructed, its location
    gets a refgetAccession mutation, and then it's identified. If the helper
    is bypassed and ``ga4gh_identify`` is called directly, the Merkle-tree
    returns the stale digest without recomputing, so the identifier doesn't
    reflect the current content. By clearing cached digests first, the helper
    ensures the identifier is always correct even if the input allele has been
    mutated since the last identification.
    """
    allele = _make_allele(start=5, end=6)
    allele.location.digest = "ga4gh:SL.stale-location-digest"
    allele.digest = "ga4gh:VA.stale-allele-digest"

    identifier = identify_allele(allele)

    assert identifier.startswith("ga4gh:VA.")
    assert identifier != "ga4gh:VA.stale-allele-digest"
    # The identifier must match a freshly-built allele with the same content,
    # not anything derived from the stale digests.
    assert identifier == identify_allele(_make_allele(start=5, end=6))


def test_normalize_and_identify_pairs_both_steps(mocker):
    """normalize_and_identify normalizes then identifies in that order."""
    allele = _make_allele()
    # Pass-through normalize so the test stays offline (no SeqRepo).
    mock_normalize = mocker.patch(
        "dcd_mapping.vrs_utils.normalize", side_effect=lambda a, **_: a
    )
    mocker.patch("dcd_mapping.vrs_utils.get_seqrepo", return_value=mocker.MagicMock())

    result = normalize_and_identify(allele)

    mock_normalize.assert_called_once()
    assert result is allele
    assert result.id is not None
    assert result.id.startswith("ga4gh:VA.")


def test_normalize_and_identify_clears_stale_digest_through_helper(mocker):
    """Stale digests survive a no-op normalize, so the identify step must
    still clear them. Verifies the pairing actually delivers the invariant.
    """
    allele = _make_allele()
    allele.location.digest = "ga4gh:SL.stale"
    allele.digest = "ga4gh:VA.stale"
    mocker.patch("dcd_mapping.vrs_utils.normalize", side_effect=lambda a, **_: a)
    mocker.patch("dcd_mapping.vrs_utils.get_seqrepo", return_value=mocker.MagicMock())

    result = normalize_and_identify(allele)

    assert result.id != "ga4gh:VA.stale"
    assert result.id.startswith("ga4gh:VA.")


def test_identify_allele_raises_when_digest_unobtainable(mocker):
    """When ``ga4gh_identify`` returns ``None`` (malformed input), surface it
    as a ValueError rather than silently stamping an unidentifiable allele.
    """
    mocker.patch("dcd_mapping.vrs_utils.ga4gh_identify", return_value=None)
    with pytest.raises(ValueError, match="Failed to compute GA4GH identifier"):
        identify_allele(_make_allele())
