"""Parse fcs files into AnnData.

Import the package::

   import readfcs

This is the complete API reference:

.. autosummary::
   :toctree: .

   read
   ReadFCS
   view
   datasets
"""

__version__ = "2.0.0"


from . import datasets
from ._core import ReadFCS, read, view
