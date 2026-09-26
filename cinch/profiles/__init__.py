"""Validated presence/type profiles and explicit legacy adapters."""

from .io import load_profile, save_profile
from .schema import LegacyMissingPolicy, StateProfile, from_legacy_alleles

__all__ = [
    "LegacyMissingPolicy",
    "StateProfile",
    "from_legacy_alleles",
    "load_profile",
    "save_profile",
]
