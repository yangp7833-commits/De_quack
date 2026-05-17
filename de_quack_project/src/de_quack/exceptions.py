class DeDuckError(Exception):
    """Base exception for all de_duck library errors."""
    pass

class ProcessingError(DeDuckError):
    """Raised when incoming dataframe columns or data types are invalid."""
    pass

class DuplicateExperimentError(DeDuckError):
    """Raised when trying to ingest an experiment_id that already exists."""
    pass

class DuplicateGeneTableError(DeDuckError):
    """Raised when the user attempts to initialize the gene table twice"""
    pass