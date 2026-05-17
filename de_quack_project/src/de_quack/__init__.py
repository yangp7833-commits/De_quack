"""de_quack package initializer."""

__version__ = "0.1.1"

from .core import de_quackling
from .viz import volcano_plot
from .exceptions import DeDuckError, ProcessingError, DuplicateExperimentError, DuplicateGeneTableError

__all__ = [
    "de_quackling",
    "volcano_plot",
    "DeDuckError",
    "ProcessingError",
    "DuplicateExperimentError",
    "DuplicateGeneTableError",
]
