"""de_quack package initializer."""

__version__ = "0.1.2"

from .core import de_quackling
from .viz import volcano_plot
from .exceptions import DeQuackError, ProcessingError, DuplicateExperimentError, DuplicateGeneTableError
from .arrow import de_arrow, de_arrows

__all__ = [
    "de_quackling",
    "volcano_plot",
    "DeQuackError",
    "ProcessingError",
    "DuplicateExperimentError",
    "DuplicateGeneTableError",
    "de_arrow",
    "de_arrows",
]
