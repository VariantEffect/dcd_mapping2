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
