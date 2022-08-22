import anndata as ad
import dateutil.parser as date_parser  # type: ignore
import fcsparser
import numpy as np
import pandas as pd
from lamin_logger import logger


def _has_numbers(inputString):
    """Check if a string contains any numbers."""
    return any(char.isdigit() for char in inputString)


def _get_spill_matrix(matrix_string: str) -> pd.DataFrame:
    """Generate a spill matrix for string.

    Code is modified from: https://github.com/whitews/FlowUtils
    Pedersen NW, Chandran PA, Qian Y, et al. Automated Analysis of Flow Cytometry
    Data to Reduce Inter-Lab Variation in the Detection of Major Histocompatibility
    Complex Multimer-Binding T Cells. Front Immunol. 2017;8:858.
    Published 2017 Jul 26. doi:10.3389/fimmu.2017.00858

    Args:
    matrix_string: str
        string value extracted from the 'spill' parameter of the FCS file

    Returns:
        Pandas.DataFrame
    """
    matrix_list = matrix_string.split(",")
    n = int(matrix_list[0])
    header = matrix_list[1 : (n + 1)]
    header = [i.strip().replace("\n", "") for i in header]
    values = [i.strip().replace("\n", "") for i in matrix_list[n + 1 :]]
    matrix = np.reshape(list(map(float, values)), (n, n))
    matrix_df = pd.DataFrame(matrix)
    matrix_df = matrix_df.rename(
        index={k: v for k, v in zip(matrix_df.columns.to_list(), header)},
        columns={k: v for k, v in zip(matrix_df.columns.to_list(), header)},
    )
    return matrix_df


class ReadFCS:
    """Read in fcs file using fcsparesr as preprocess the metadata.

    Args:
        filepath: str or Path
            location of fcs file to parse
        comp_matrix: str
            csv file containing compensation matrix (optional, not required if a
            spillover matrix is already linked to the file)
    """

    def __init__(self, filepath):
        # Added to allow reading in Path
        self._meta_raw, self._data = fcsparser.parse(filepath, reformat_meta=True)

        self._meta = {
            (k.replace("$", "") if _has_numbers(k) else k.replace("$", "").lower()): v
            for k, v in self._meta_raw.items()
        }
        self.spill_txt = None
        self.spill = None

        # header
        self._meta["header"] = self.meta["__header__"]
        self._meta["header"]["FCS format"] = self.meta["__header__"][
            "FCS format"
        ].decode()

        # channels
        self._meta["channels"] = self.meta["_channels_"]

        # date
        try:
            self._meta["processing_date"] = date_parser.parse(
                self.meta["date"] + " " + self.meta["etim"]
            ).isoformat()
        except (KeyError, date_parser.ParserError):
            self._meta["processing_date"] = None

        # compensation matrix
        spill_list = [
            self.meta.get(key)
            for key in ["spill", "spillover"]
            if self.meta.get(key) is not None
        ]
        if len(spill_list) > 0:
            self.spill_txt = spill_list[0]
            if len(self.spill_txt) > 0:
                self._meta["spill"] = _get_spill_matrix(self.spill_txt)
        else:
            logger.warning(
                "No spillover matrix found, please provide path to relevant csv file with 'comp_matrix' argument if compensation is necessary"  # noqa
            )

    @property
    def meta(self):
        """Metadata."""
        return self._meta

    @property
    def data(self):
        """Data matrix."""
        return self._data

    def compensate(self):
        """Apply compensation to event data."""
        assert (
            self.meta["spill"] is not None
        ), f"Unable to locate spillover matrix, please provide a compensation matrix"  # noqa
        channel_idx = [
            i
            for i, (idx, row) in enumerate(self._meta["channels"].iterrows())
            if row["$PnS"] not in ["", " "]
        ]

        channel_idx = [
            i
            for i, (idx, row) in enumerate(self._meta["channels"].iterrows())
            if all([z not in row["$PnN"].lower() for z in ["fsc", "ssc", "time"]])
            and row["$PnN"] in self.meta["spill"].columns  # noqa
        ]

        comp_data = self.data.iloc[:, channel_idx]
        comp_data = np.linalg.solve(self.meta["spill"].values.T, comp_data.T).T
        self._data.iloc[:, channel_idx] = comp_data

    def to_anndata(self, reindex=True) -> ad.AnnData:
        """Convert the FCSFile instance to an AnnData.

        Args:
            reindex: variables will be reindexed with marker names if possible otherwise
                channels
        Returns:
            an AnnData object
        """
        channels_mapping = {
            "$PnN": "channel",
            "$PnS": "marker",
        }

        adata = ad.AnnData(
            self.data,
            var=pd.DataFrame(
                self._meta["channels"].rename(columns=channels_mapping)
            ).set_index("channel"),
        )

        if reindex:
            adata.var = adata.var.reset_index()
            adata.var.index = np.where(
                adata.var["marker"].isin(["", " "]),
                adata.var["channel"],
                adata.var["marker"],
            )
            mapper = pd.Series(adata.var.index, index=adata.var["channel"])
            if self.meta["spill"] is not None:
                self._meta["spill"].index = self.meta["spill"].index.map(mapper)
                self._meta["spill"].columns = self.meta["spill"].columns.map(mapper)

        adata.uns["meta"] = self.meta

        return adata


def read(filepath, reindex=True) -> ad.AnnData:
    """Read in fcs file as AnnData.

    Args:
        filepath: str or Path
            location of fcs file to parse
        reindex: bool
            variables will be reindexed with marker names if possible otherwise
            channels
    Returns:
        an AnnData object
    """
    fcsfile = ReadFCS(filepath)
    return fcsfile.to_anndata(reindex=reindex)
