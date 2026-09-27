"""Unified statistical APIs with explicit scientific meanings.

The frozen v1 implementation remains available under :mod:`cinch.frozen_v1`.
This package does not alias weighted and unweighted estimators.
"""

from .epidis import directional_epidis
from .multiple_testing import bh_adjust
from .weighted_mi import weighted_binary_information

__all__ = ["bh_adjust", "directional_epidis", "weighted_binary_information"]
