"""gotools: fast conserved-site and variant detection from MAF alignments.

Three command-line tools, ``go2fix``, ``go2var``, and ``addpro``, are
installed as console scripts when the package is installed. The same
functionality is available as a Python API; the supported entry points
are re-exported below.
"""

from gotools.addpro import addpro_gtf_to_bed
from gotools.go2fix import go2fix_optimized, go2fix_single
from gotools.go2var import go2var_sorted_optimized

__version__ = "0.2.0"

__all__ = [
    "__version__",
    "addpro_gtf_to_bed",
    "go2fix_optimized",
    "go2fix_single",
    "go2var_sorted_optimized",
]
