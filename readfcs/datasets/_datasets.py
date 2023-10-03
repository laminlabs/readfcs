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


def Oetjen18_t1() -> str:
    """A fcs file with T cell panel from Oetjen18.

    Original file name: 2-13-17 T cell Panel_T_E_G05_004.fcs

    Downloaded from: http://flowrepository.org/id/FR-FCM-ZYQ9

    Reference: https://insight.jci.org/articles/view/124928

    """
    url = "https://lamindb-dev-datasets.s3.amazonaws.com/.lamindb/DBNEczSgBui0bbzBXMGH.fcs"  # noqa
    path_data, _ = urlretrieve(url, "oetjen18_t1.fcs")

    return path_data


def Oetjen18_t2() -> str:
    """A fcs file with T cell panel from Oetjen18.

    Original file name: 2-13-17 T cell Panel_T_Q_G04_003.fcs

    Downloaded from: http://flowrepository.org/id/FR-FCM-ZYQ9

    Reference: https://insight.jci.org/articles/view/124928
    """
    url = "https://lamindb-dev-datasets.s3.amazonaws.com/.lamindb/ckbFvcuxG4tln7OIg3ml.fcs"  # noqa
    path_data, _ = urlretrieve(url, "oetjen18_t2.fcs")

    return path_data


def Oetjen18_dc() -> str:
    """A fcs file with DC panel from Oetjen18.

    Original file name: DC Panel_DC_Q_F04_003.fcs

    Downloaded from: http://flowrepository.org/id/FR-FCM-ZYQ9

    Reference: https://insight.jci.org/articles/view/124928
    """
    url = "https://lamindb-dev-datasets.s3.amazonaws.com/.lamindb/yYLIWRT3sg4E3NUFigpB.fcs"  # noqa
    path_data, _ = urlretrieve(url, "oetjen18_dc.fcs")

    return path_data
