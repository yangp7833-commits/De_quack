"""de_quack package initializer."""

__version__ = "0.1.2"

from .arrow import DeArrow, DeArrows
from .core import DeQuackling
from .exceptions import (
    DeQuackError,
    DuplicateExperimentError,
    DuplicateGeneTableError,
    ProcessingError,
)
from .viz import volcano_plot

__all__ = [
    "DeArrow",
    "DeArrows",
    "DeQuackError",
    "DeQuackling",
    "DuplicateExperimentError",
    "DuplicateGeneTableError",
    "ProcessingError",
    "volcano_plot",
]
