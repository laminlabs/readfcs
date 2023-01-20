import re
from pathlib import Path
from typing import Dict, Union

import anndata as ad
import fcsparser
import numpy as np
import pandas as pd


def _channels_to_df(meta_raw: dict) -> dict:
    """Format channels into a DataFrame.

    Args:
    meta_raw: dict
        original metadata

    Returns:
        a dict of metadata with "channels"
    """
    meta = {}
    channel_groups: Dict = {}

    # channel groups are $PnB, $PnS, $PnN...
    for k, v in meta_raw.items():
        # Get all fields with $PnX pattern
        if re.match("^\$P\d+[A-Z]$", k):  # noqa
            group_key = f"$Pn{k[-1]}"
            if group_key not in channel_groups:
                channel_groups[group_key] = []
            # The numeric index n
            idx = int(k.lstrip("$P")[:-1])
            channel_groups[group_key].append((idx, v))
        else:
            meta[k] = v

    # format channels into a dataframe
    k = "$PnN"
    df_groups = pd.DataFrame(channel_groups.get(k), columns=["n", k]).set_index("n")

    for k, group in channel_groups.items():
        if k == "$PnN":
            continue
        df_group = pd.DataFrame(channel_groups.get(k), columns=["n", k]).set_index("n")
        df_groups = df_groups.join(df_group)

    # reorder df_groups columns
    df_groups.insert(0, "$PnN", df_groups.pop("$PnN"))
    if "$PnS" in df_groups.columns:
        df_groups.insert(1, "$PnS", df_groups.pop("$PnS"))
    meta["channels"] = df_groups

    return meta


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
        data_set: int
            Index of retrieved data set in the fcs file.
    """

    def __init__(self, filepath: Union[str, Path], data_set: int = 0) -> None:
        # No metadata formating using fcsparser
        self._meta_raw, self._data = fcsparser.parse(
            filepath, data_set=data_set, channel_naming="$PnN"
        )

        # Format channels into a dataframe as `self.meta["channels"]`
        self._meta = _channels_to_df(self._meta_raw)

        # header
        self._meta["header"] = self.meta["__header__"]
        self._meta["header"]["FCS format"] = self.meta["__header__"][
            "FCS format"
        ].decode()

        # compensation matrix
        self.spill_txt = None
        self.spill = None
        if "SPILL" in self.meta:
            self._meta["spill"] = self._meta.pop("SPILL")
        if "SPILLOVER" in self.meta:
            self._meta["spillover"] = self._meta.pop("SPILLOVER")

        spill_list = [
            self.meta.get(key)
            for key in ["spill", "spillover"]
            if self.meta.get(key) is not None
        ]
        if len(spill_list) > 0:
            self.spill_txt = spill_list[0]
            if len(self.spill_txt) > 0:  # type:ignore
                self._meta["spill"] = _get_spill_matrix(self.spill_txt)  # type:ignore

    @property
    def header(self) -> dict:
        """Header."""
        return self._meta["header"]

    @property
    def meta(self) -> dict:
        """Metadata."""
        return self._meta

    @property
    def data(self) -> pd.DataFrame:
        """Data matrix."""
        return self._data

    def compensate(self) -> None:
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
        self._data[self._data.columns[channel_idx]] = comp_data

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
        if any([i for i in ["$PnN", "$PnS"] if i not in self.meta["channels"].columns]):
            raise AssertionError(
                "$PnN or $PnS field not found in the file!\nPlease check your file"
                " content with `readfcs.view`!"
            )

        # AnnData only allows str index
        var = self._meta["channels"]
        X = self._data
        var.index = var.index.astype(str)
        X.columns = var.index
        X.index = X.index.astype(str)

        # convert list columns to str so it saves properly
        for k in var.columns:
            if isinstance(var[k].iloc[0], list):
                var[k] = var[k].map(repr)

        # create anndata with channel indexing variables
        adata = ad.AnnData(
            X,
            var=var.rename(columns=channels_mapping),
        )

        # by default, we index variables with marker
        # use channels for non-marker channels
        if reindex:
            adata.var = adata.var.reset_index()
            adata.var.index = np.where(
                adata.var["marker"].isin(["", " ", np.nan]),
                adata.var["channel"],
                adata.var["marker"],
            )
            mapper = pd.Series(adata.var.index, index=adata.var["channel"])
            if self.meta.get("spill") is not None:
                self._meta["spill"].index = self.meta["spill"].index.map(mapper)
                self._meta["spill"].columns = self.meta["spill"].columns.map(mapper)

        # write metadata into adata.uns
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


def view(filepath: Union[str, Path], data_set: int = 0):
    """Read in file content without preprocessing for debugging.

    Args:
        filepath: str or Path
            location of fcs file to parse
        data_set: int
            Index of retrieved data set in the fcs file.

    Returns:
        a tuple of (data, metadata)
        - data is a DataFrame
        - metadata is a dictionary

    See `fcsparser.parse`: https://github.com/eyurtsev/fcsparser
    """
    meta, data = fcsparser.parse(filepath, data_set=data_set, channel_naming="$PnN")
    return meta, data
