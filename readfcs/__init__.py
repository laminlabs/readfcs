"""Parse fcs files into AnnData.

Import the package::

   import readfcs

This is the complete API reference:

.. autosummary::
   :toctree: .

   read
   FCSFile
   datasets
"""

__version__ = "0.1.4"

from . import datasets
from ._core import FCSFile, read
