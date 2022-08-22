from pathlib import PosixPath

import anndata as ad
import dateutil.parser as date_parser  # type: ignore
import flowio
import numpy as np
import pandas as pd
from lamin_logger import logger


def chunks(df_list: list, n: int) -> pd.DataFrame:
    """Yield successive n-sized chunks from l.

    ref: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks # noqa

    Args:
        df_list: list
            list of DataFrames to generated 'chunks' from
        n: int
            number of chunks to generate

    Returns:
        Yields successive n-sized DataFrames
    """
    for i in range(0, len(df_list), n):
        yield df_list[i : i + n]


def _get_spill_matrix(matrix_string: str) -> pd.DataFrame:
    """Generate dataframe for spillover matrix.

    Generate pandas dataframe for the fluorochrome spillover matrix used for
    compensation calc

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


def _get_channel_mappings(fluoro_dict: dict) -> list:
    """Generate a list of dicts for fluorochrome mappings.

    Code is modified from cytopy:
    https://github.com/burtonrj/CytoPy/blob/master/cytopy/data/read_write.py

    Args:
        fluoro_dict: dict
            dictionary object from the channels param of the fcs file

    Returns:
        List of dict obj with keys 'channel' and 'marker'.
    """
    fm = [(int(k), x) for k, x in fluoro_dict.items()]
    fm = [x[1] for x in sorted(fm, key=lambda x: x[0])]
    mappings = []
    for fm_ in fm:
        channel = fm_["PnN"].replace("_", "-")  # type: ignore
        marker = ""
        if "PnS" in fm_.keys():  # type: ignore
            marker = fm_["PnS"].replace("_", "-")  # type: ignore
        mappings.append({"channel": channel, "marker": marker})
    return mappings


class FCSFile:
    """Utilising FlowIO to generate an object for representing an FCS file.

    Code is modified from cytopy
    https://github.com/burtonrj/CytoPy/blob/master/cytopy/data/read_write.py
    to:
    - allow reading in Path object
    - catch ValueError when reading in fcs using flowio.FlowData
    - catch ParserError in processing date
    - fix the shape mismatch in .compensate()

    Args:
        filepath: str or Path
            location of fcs file to parse
        comp_matrix: str
            csv file containing compensation matrix (optional, not required if a
            spillover matrix is already linked to the file)
    """

    def __init__(self, filepath, comp_matrix=None):
        # Added to allow reading in Path
        if isinstance(filepath, PosixPath):
            filepath = filepath.as_posix()
        fcs = flowio.FlowData(filepath)
        self._fcs = fcs

        # metadata from the fcs.text
        self.filename = fcs.text.get("fil", "Unknown_filename")
        self.sys = fcs.text.get("sys", "Unknown_system")
        self.total_events = int(fcs.text.get("tot", 0))
        self.tube_name = fcs.text.get("tube name", "Unknown")
        self.exp_name = fcs.text.get("experiment name", "Unknown")
        self.cytometer = fcs.text.get("cyt", "Unknown")
        self.creator = fcs.text.get("creator", "Unknown")
        self.operator = fcs.text.get("export user name", "Unknown")
        self.cst_pass = fcs.text.get("cst setup status", "FAIL")
        self.threshold = "Unknown"
        self.spill_txt = None
        self.spill = None

        # additional metadata
        self.header = fcs.header
        self.channel_mappings = _get_channel_mappings(fcs.channels)

        # data
        self.events = fcs.events
        self.data = np.reshape(
            np.array(fcs.events, dtype=np.float32), (-1, fcs.channel_count)
        )

        # channels
        if "threshold" in fcs.text.keys():
            self.threshold = [
                {"channel": c, "threshold": v}
                for c, v in chunks(fcs.text["threshold"].split(","), 2)
            ]

        try:
            self.processing_date = date_parser.parse(
                fcs.text["date"] + " " + fcs.text["etim"]
            ).isoformat()
        except (KeyError, date_parser.ParserError):
            self.processing_date = "Unknown"

        # compensation matrix
        if comp_matrix is not None:
            self.spill = pd.read_csv(comp_matrix)
        else:
            spill_list = [
                fcs.text.get(key)
                for key in ["spill", "spillover"]
                if fcs.text.get(key) is not None
            ]
            if len(spill_list) > 0:
                self.spill_txt = spill_list[0]
                if len(self.spill_txt) > 0:
                    self.spill = _get_spill_matrix(self.spill_txt)
            else:
                logger.warning(
                    "No spillover matrix found, please provide path to relevant csv file with 'comp_matrix' argument if compensation is necessary"  # noqa
                )

    def compensate(self):
        """Apply compensation to event data."""
        assert (
            self.spill is not None
        ), f"Unable to locate spillover matrix, please provide a compensation matrix"  # noqa
        channel_idx = [
            i for i, x in enumerate(self.channel_mappings) if x["marker"] != ""
        ]

        channel_idx = [
            i
            for i, x in enumerate(self.channel_mappings)
            if all([z not in x["channel"].lower() for z in ["fsc", "ssc", "time"]])
            and x["channel"] in self.spill.columns  # noqa
        ]

        comp_data = self.data[:, channel_idx]
        comp_data = np.linalg.solve(self.spill.values.T, comp_data.T).T
        self.data[:, channel_idx] = comp_data

    def write_fcs(self, filename, metadata=None):
        """Export FCSFile instance as a new FCS file.

        Args:
            filename: name of exported FCS file
            metadata: an optional dictionary for adding metadata keywords/values

        """
        self._fcs.write_fcs(filename=filename, metadata=metadata)

    def to_anndata(self, reindex=True) -> ad.AnnData:
        """Convert the FCSFile instance to an AnnData.

        Args:
            reindex: variables will be reindexed with marker names if possible otherwise
                channels
        Returns:
            an AnnData object
        """
        adata = ad.AnnData(
            self.data,
            var=pd.DataFrame(self.channel_mappings).set_index("channel"),
        )
        if reindex:
            adata.var = adata.var.reset_index()
            adata.var.index = np.where(
                adata.var["marker"].isin(["", " "]),
                adata.var["channel"],
                adata.var["marker"],
            )
            mapper = pd.Series(adata.var.index, index=adata.var["channel"])
            if self.spill is not None:
                self.spill.index = self.spill.index.map(mapper)
                self.spill.columns = self.spill.columns.map(mapper)
        meta = {
            "filename": self.filename,
            "sys": self.sys,
            "header": self.header,
            "total_events": self.total_events,
            "tube_name": self.tube_name,
            "exp_name": self.exp_name,
            "cytometer": self.cytometer,
            "creator": self.creator,
            "operator": self.operator,
            "cst_pass": self.cst_pass,
            "threshold": self.threshold[0],
            "processing_date": self.processing_date,
            "spill": self.spill,
        }
        adata.uns["meta"] = meta

        return adata


def read(filepath, comp_matrix=None, reindex=True) -> ad.AnnData:
    """Read in fcs file as AnnData.

    Args:
        filepath: str or Path
            location of fcs file to parse
        comp_matrix: str
            csv file containing compensation matrix (optional, not required if a
            spillover matrix is already linked to the file)
        reindex: bool
            variables will be reindexed with marker names if possible otherwise
            channels
    Returns:
        an AnnData object
    """
    fcsfile = FCSFile(filepath, comp_matrix=comp_matrix)
    return fcsfile.to_anndata(reindex=reindex)
