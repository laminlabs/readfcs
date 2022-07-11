"""Parse fcs files into AnnData.

Import the package::

   import fcsreader

This is the complete API reference:

.. autosummary::
   :toctree: .

   FCSFile
"""

__version__ = "0.1.0"

from . import datasets
from ._core import FCSFile
