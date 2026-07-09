"""Exceptions for DCD Mapping Module"""


### Custom Exceptions for VRS Mapping Errors ###


class VrsMapError(Exception):
    """Raise in case of generic VRS mapping errors."""


class UnsupportedReferenceSequenceNameSpaceError(ValueError):
    """Raised when a reference sequence name space is not supported."""


class MissingSequenceIdError(ValueError):
    """Raised when a sequence ID is not provided."""


class UnsupportedReferenceSequencePrefixError(ValueError):
    """Raised when a reference sequence prefix is not supported."""


### Custom Exceptions for Alignment Errors ###


class AlignmentError(ValueError):
    """Raise when errors encountered during alignment."""


class BlatNotFoundError(AlignmentError):
    """Raise when BLAT binary appears to be missing."""


### Custom Exceptions for Data Lookup Errors ###


class DataLookupError(Exception):
    """Raise for misc. issues related to resource acquisition/lookup."""


### Custom Exceptions for MaveDB Data Errors ###


class ScoresetNotSupportedError(ValueError):
    """Raise when a score set cannot be mapped because it has characteristics that are not currently supported."""


### Custom Exceptions for Resource Acquisition Errors ###


class ResourceAcquisitionError(ValueError):
    """Raise when resource acquisition fails."""


### Custom Exceptions for Transcript Selection Errors ###


class TxSelectError(ValueError):
    """Raise for transcript selection failure."""


class NoCodingTranscriptError(TxSelectError):
    """Raise when a protein-coding target has no resolvable coding transcript.

    Distinct from the regulatory/non-coding case (which returns ``None`` from
    transcript selection because no coding transcript is expected): this signals
    a coding target for which selection *should* have produced a transcript but
    could not -- no resolvable gene symbol, no MANE/compatible transcript for the
    gene, or a projection that does not land cleanly on the selected transcript.
    Downstream records this as a recoverable skip, distinct from "no protein
    consequence exists."
    """
