[![Stars](https://img.shields.io/github/stars/laminlabs/readfcs?logo=GitHub&color=yellow)](https://github.com/laminlabs/readfcs)
[![codecov](https://codecov.io/gh/laminlabs/readfcs/branch/main/graph/badge.svg?token=6A5PYRX809)](https://codecov.io/gh/laminlabs/readfcs)
[![pypi](https://img.shields.io/pypi/v/readfcs?color=blue&label=pypi%20package)](https://pypi.org/project/readfcs)
[![doi](https://img.shields.io/badge/doi-10.56528%2Frfcs-lightgrey)](https://doi.org/10.56528/rfcs)

# readfcs: Read FCS files

Lightweight Python library to load data and metadata from Flow Cytometry Standard (FCS) files into `DataFrame` and `AnnData` objects.

<br>

Install: ![pyversions](https://img.shields.io/pypi/pyversions/readfcs)

```bash
$ pip install readfcs
```

Get started:

```python
import readfcs

readfcs.read()  # Read FCS files into `AnnData`

readfcs.ReadFCS()  # Debug possible integrity issues in your FCS files
```

Read the [docs](https://readfcs.lamin.ai/) or read the [blog post](https://lamin.ai/blog/2022/readfcs).
