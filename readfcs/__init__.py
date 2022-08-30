"""Parse fcs files into AnnData.

Import the package::

   import readfcs

This is the complete API reference:

.. autosummary::
   :toctree: .

   read
   ReadFCS
   datasets
"""

__version__ = "1.0.1"

from . import datasets
from ._core import ReadFCS, read
