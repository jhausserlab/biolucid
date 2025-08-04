"""
bioLUCID - A tool for batch effect analysis in single-cell RNA sequencing data.
"""


from . import core
from . import visualization
from . import utils
from . import config

__name__ = "biolucid"

__version__ = "1.1.0"

__all__ = [
    "core",
    "visualization",
    "utils",
    "config",
] 