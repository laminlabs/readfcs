from urllib.request import urlretrieve


def example() -> str:
    """An example fcs file.

    With a spill matrix.

    Originally from:
    https://github.com/whitews/FlowIO/blob/master/examples/fcs_files/100715.fcs
    """
    url = "https://lamindb-test.s3.amazonaws.com/example.fcs"
    path_data, _ = urlretrieve(url, "example.fcs")

    return path_data


def flowio_3FITC_4PE_004() -> str:
    """An example fcs file.

    Originally from:
    https://github.com/whitews/FlowIO/blob/master/examples/fcs_files/3FITC_4PE_004.fcs
    """
    url = "https://github.com/whitews/FlowIO/blob/master/examples/fcs_files/3FITC_4PE_004.fcs"  # noqa
    path_data, _ = urlretrieve(url, "3FITC_4PE_004.fcs")

    return path_data
