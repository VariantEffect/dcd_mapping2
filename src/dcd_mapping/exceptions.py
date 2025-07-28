"""Exceptions for DCD Mapping Module"""


class VrsMapError(Exception):
    """Raise in case of generic VRS mapping errors."""


class UnsupportedReferenceSequenceNameSpaceError(ValueError):
    """Raised when a reference sequence name space is not supported."""


class MissingSequenceIdError(ValueError):
    """Raised when a sequence ID is not provided."""


class UnsupportedReferenceSequencePrefixError(ValueError):
    """Raised when a reference sequence prefix is not supported."""
