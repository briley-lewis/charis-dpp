"""CHARIS Data Processing Pipeline (Python port).

Sub-packages mirror the IDL layout: ``prelims``, ``skysub``, ``reg``,
``specphotcal``, ``filt``, ``psfsub``, ``partialsub``, ``planetspec``,
``pdi``, ``tools``.  Shared infrastructure lives in ``config``, ``io``,
``geometry``, ``image`` and ``astro``.
"""

import warnings
from importlib.metadata import PackageNotFoundError, version

warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="Mean of empty slice")

try:
    __version__ = version("charis-dpp")
except PackageNotFoundError:  # pragma: no cover - source checkout
    __version__ = "0.0.0+src"

__all__ = ["__version__"]
