"""VRS allele identification helpers.

Centralizes the digest-correctness invariant for GA4GH VRS alleles: the
``ga4gh_identify`` Merkle tree caches sub-object digests on the object after
first identification, so any subsequent mutation (refgetAccession swap,
normalization, state coercion) leaves a stale id unless the cached digests
are cleared first. All allele identification in dcd_mapping must route
through :func:`identify_allele` so the digest is always recomputed from
current content.
"""

from ga4gh.core import ga4gh_identify
from ga4gh.vrs._internal.models import Allele, SequenceLocation
from ga4gh.vrs.normalize import normalize

from dcd_mapping.lookup import get_seqrepo


def identify_allele(allele: Allele) -> str:
    """Clear cached digests and return a fresh GA4GH identifier for *allele*.

    ``ga4gh_identify`` is a Merkle tree: it calls ``get_or_create_digest`` on
    sub-objects, returning any cached value without recomputing. Clearing both
    the location digest and the allele digest first ensures the id is derived
    from current content, not from a value set before a refgetAccession
    mutation or normalization.

    Clearing the digests alone isn't enough: ``in_place="default"`` only fills
    an *empty* ``id``, so an allele that already has one is handed that value
    straight back. Every allele reaching this helper already has one —
    ``AlleleTranslator`` defaults to ``identify=True`` and is called with
    ``do_normalize=False``, so the id is minted over *un-normalized* content.
    ``in_place="always"`` forces recomputation from the freshly cleared digest.

    This only matters where normalization moves the span, which is why it went
    unnoticed: a substitution normalizes to itself, so its pre-normalization id
    stays correct by coincidence. Deletions and duplications in repeat regions
    expand leftward, and those are the ones that end up mislabelled.

    The location must be identified separately: ``ga4gh_identify`` writes only
    the ``id`` of the object it's handed, and for sub-objects calls
    ``get_or_create_digest``, never ``get_or_create_ga4gh_identifier``. Clearing
    the location's digest fixes what the *allele* id derives from, but leaves
    ``location.id`` holding the pre-normalization value unless it too is
    explicitly recomputed.
    """
    if isinstance(allele.location, SequenceLocation):
        allele.location.digest = None
        allele.location.id = None
        ga4gh_identify(allele.location, in_place="always")

    allele.digest = None
    digest = ga4gh_identify(allele, in_place="always")
    if digest is None:
        raise ValueError("Failed to compute GA4GH identifier for allele")  # noqa: EM101

    return digest


def normalize_and_identify(allele: Allele) -> Allele:
    """Normalize *allele* and stamp it with a freshly computed GA4GH digest.

    Pairs the two finalize steps every VRS allele construction path needs.
    Routing identification through :func:`identify_allele` (rather than
    ``ga4gh_identify`` directly) is the invariant that protects against the
    Merkle-tree's stale-digest behavior after mutation -- so any allele
    construction site that bypasses this helper risks reintroducing the
    stale-digest bug.
    """
    allele = normalize(allele, data_proxy=get_seqrepo())
    allele.id = identify_allele(allele)
    return allele
