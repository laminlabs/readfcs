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

__version__ = "1.0.4"

# prints warning of python versions
from lamin_logger import py_version_warning

py_version_warning("3.7", "3.10")

from . import datasets
from ._core import ReadFCS, read
