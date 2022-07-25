"""Parse fcs files into AnnData.

Import the package::

   import readfcs

This is the complete API reference:

.. autosummary::
   :toctree: .

   FCSFile
"""

__version__ = "0.1.1"

from . import datasets
from ._core import FCSFile
