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

__version__ = "1.1.8"


from . import datasets
from ._core import ReadFCS, read, view
