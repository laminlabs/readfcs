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

__version__ = "0.1.7"

from . import datasets
from ._core import ReadFCS, read
