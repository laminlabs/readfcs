import json
import os
from multiprocessing import Pool, cpu_count
from pathlib import PosixPath
from typing import Dict, List, Optional

import anndata as ad
import dateutil.parser as date_parser  # type: ignore
import flowio
import numpy as np
import pandas as pd


def filter_fcs_files(
    fcs_dir: str, exclude_comps: bool = True, exclude_dir: str = "DUPLICATES"
) -> list:
    """Dict -> fcs file paths.

    Given a directory, return file paths for all fcs files in directory and
    subdirectories contained within.

    Args:
        fcs_dir: str
            path to directory for search
        exclude_comps: bool
            if True, compensation files will be ignored (note: function
            searches for 'comp' in file name for exclusion)
        exclude_dir: str (default = 'DUPLICATES')
            Will ignore any directories with this name

    Returns:
        list of fcs file paths
    """
    fcs_files: List[str] = []
    for root, _, files in os.walk(fcs_dir):
        if os.path.basename(root) == exclude_dir:
            continue
        if exclude_comps:
            fcs = [
                f
                for f in files
                if f.lower().endswith(".fcs") and f.lower().find("comp") == -1
            ]
        else:
            fcs = [f for f in files if f.lower().endswith(".fcs")]
        fcs = [os.path.join(root, f) for f in fcs]
        fcs_files = fcs_files + fcs
    return fcs_files


def get_fcs_file_paths(
    fcs_dir: str,
    control_names: list,
    ctrl_id: str,
    ignore_comp: bool = True,
    exclude_dir: str = "DUPLICATE",
) -> dict:
    """Generate a standard dictionary object of fcs files in given directory.

    Args:
        fcs_dir: str
            target directory for search
        control_names: list
            names of expected control files (names must appear in filenames)
        ctrl_id: str
            global identifier for control file e.g. 'FMO' (must appear in filenames)
        ignore_comp: bool, (default=True)
            If True, files with 'compensation' in their name will be ignored
            (default = True)
        exclude_dir: str (default = 'DUPLICATES')
            Will ignore any directories with this name

    Returns:
        standard dictionary of fcs files contained in target directory
    """
    file_tree: Dict = dict(primary=[], controls={})
    fcs_files = filter_fcs_files(
        fcs_dir, exclude_comps=ignore_comp, exclude_dir=exclude_dir
    )
    ctrl_files = [f for f in fcs_files if f.find(ctrl_id) != -1]
    primary = [f for f in fcs_files if f.find(ctrl_id) == -1]
    for c_name in control_names:
        matched_controls = list(filter(lambda x: x.find(c_name) != -1, ctrl_files))
        if not matched_controls:
            print(f"Warning: no file found for {c_name} control")
            continue
        if len(matched_controls) > 1:
            print(f"Warning: multiple files found for {c_name} control")
        file_tree["controls"][c_name] = matched_controls

    if len(primary) > 1:
        print(
            "Warning! Multiple non-control (primary) files found in directory. Check"
            " before proceeding."
        )
    file_tree["primary"] = primary
    return file_tree


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


def fcs_mappings(path: str) -> Optional[list]:
    """Fetch channel mappings from fcs file.

    Args:
        path: str
            path to fcs file

    Returns:
        List of channel mappings. Will return None if file fails to load.
    """
    try:
        fo = FCSFile(path)
    except ValueError as e:
        print(f"Failed to load file {path}; {e}")
        return None
    return fo.channel_mappings


def explore_channel_mappings(fcs_dir: str, exclude_comps: bool = True) -> list:
    """Dict -> channel/marker mappings.

    Given a directory, explore all fcs files and find all permutations of
    channel/marker mappings.

    Args:
    fcs_dir: str
        root directory to search
    exclude_comps: bool, (default=True)
        exclude compentation files (must have 'comp' in filename)

    Returns:
        list of all unique channel/marker mappings
    """
    fcs_files = filter_fcs_files(fcs_dir, exclude_comps)
    with Pool(cpu_count()) as pool:
        mappings = list(pool.map(fcs_mappings, fcs_files))  # type: ignore
        mappings = list(pool.map(json.dumps, mappings))  # type: ignore
    return [json.loads(x) for x in mappings]  # type: ignore


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

    Generates a list of dictionary objects that describe the fluorochrome
    mappings in this FCS file.

    Args:
        fluoro_dict: dict
            dictionary object from the channels param of the fcs file

    Returns:
        List of dict obj with keys 'channel' and 'marker'. Use to map fluorochrome channels to corresponding marker # noqa
    """
    fm = [(int(k), x) for k, x in fluoro_dict.items()]
    fm = [x[1] for x in sorted(fm, key=lambda x: x[0])]
    mappings = []
    for fm_ in fm:
        channel = fm_["PnN"].replace("_", "-")  # type: ignore
        if "PnS" in fm_.keys():  # type: ignore
            marker = fm_["PnS"].replace("_", "-")  # type: ignore
        else:
            marker = ""
        mappings.append({"channel": channel, "marker": marker})
    return mappings


class FCSFile:
    """Utilising FlowIO to generate an object for representing an FCS file.

    Code is modified from cytopy data.read_write.py to:
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
        try:
            fcs = flowio.FlowData(filepath)
        except ValueError:
            fcs = flowio.FlowData(filepath, ignore_offset_error=True)
        self._fcs = fcs
        self.filename = fcs.text.get("fil", "Unknown_filename")
        self.sys = fcs.text.get("sys", "Unknown_system")
        self.total_events = int(fcs.text.get("tot", 0))
        self.tube_name = fcs.text.get("tube name", "Unknown")
        self.exp_name = fcs.text.get("experiment name", "Unknown")
        self.cytometer = fcs.text.get("cyt", "Unknown")
        self.creator = fcs.text.get("creator", "Unknown")
        self.operator = fcs.text.get("export user name", "Unknown")
        self.channel_mappings = _get_channel_mappings(fcs.channels)
        self.cst_pass = False
        self.data = fcs.events
        self.event_data = np.reshape(
            np.array(fcs.events, dtype=np.float32), (-1, fcs.channel_count)
        )
        if "threshold" in fcs.text.keys():
            self.threshold = [
                {"channel": c, "threshold": v}
                for c, v in chunks(fcs.text["threshold"].split(","), 2)
            ]
        else:
            self.threshold = "Unknown"
        try:
            self.processing_date = date_parser.parse(
                fcs.text["date"] + " " + fcs.text["etim"]
            ).isoformat()
        except (KeyError, date_parser.ParserError):
            self.processing_date = "Unknown"
        if comp_matrix is not None:
            self.spill = pd.read_csv(comp_matrix)
            self.spill_txt = None
        else:
            if "spill" in fcs.text.keys():
                self.spill_txt = fcs.text["spill"]

            elif "spillover" in fcs.text.keys():
                self.spill_txt = fcs.text["spillover"]
            else:
                self.spill_txt = None
            if self.spill_txt is not None:
                if (len(self.spill_txt)) < 1:
                    print(
                        """Warning: no spillover matrix found, please provide
                    path to relevant csv file with 'comp_matrix' argument if compensation is necessary"""  # noqa
                    )
                    self.spill = None
                else:
                    self.spill = _get_spill_matrix(self.spill_txt)
            else:
                self.spill = None
        if "cst_setup_status" in fcs.text:
            if fcs.text["cst setup status"] == "SUCCESS":
                self.cst_pass = True

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

        comp_data = self.event_data[:, channel_idx]
        comp_data = np.linalg.solve(self.spill.values.T, comp_data.T).T
        self.event_data[:, channel_idx] = comp_data

    def write_fcs(self, filename, metadata=None):
        """Export FCSFile instance as a new FCS file.

        Args:
            filename: name of exported FCS file
            metadata: an optional dictionary for adding metadata keywords/values

        """
        self._fcs.write_fcs(filename=filename, metadata=metadata)

    def to_anndata(self, reindex=False) -> ad.AnnData:
        """Convert the FCSFile instance to an AnnData.

        Args:
            reindex: variables will be reindexed with marker names if possible otherwise
                channels
        Returns:
            an AnnData object
        """
        adata = ad.AnnData(
            self.event_data,
            var=pd.DataFrame(self.channel_mappings).set_index("channel"),
        )
        if reindex:
            adata.var = adata.var.reset_index()
            adata.var.index = np.where(
                adata.var["marker"] == " ", adata.var["channel"], adata.var["marker"]
            )
        meta = {
            "filename": self.filename,
            "sys": self.sys,
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


def read(filepath, comp_matrix=None, reindex=False) -> ad.AnnData:
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
