"""gotools: fast conserved-site and variant detection from MAF alignments.

Two command-line tools, ``go2fix`` and ``go2var``, are installed as console
scripts when the package is installed. The same functionality is available as
a Python API; the supported entry points are re-exported below.
"""

from gotools.go2fix import go2fix_optimized, go2fix_single
from gotools.go2var import go2var_sorted_optimized

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "go2fix_optimized",
    "go2fix_single",
    "go2var_sorted_optimized",
]
