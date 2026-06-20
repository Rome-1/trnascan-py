"""trnascan-py: a Python interface to tRNAscan-SE.

v1 binds Infernal (`cmsearch`) as the search engine and establishes a
differential oracle against reference tRNAscan-SE 2.0. v2 will accelerate the
covariance-model search core natively (JAX/GPU).
"""

from trnascan_py.models import Strand, TRNAHit
from trnascan_py.pipeline import scan

__all__ = ["Strand", "TRNAHit", "scan", "__version__"]
__version__ = "0.1.0"
