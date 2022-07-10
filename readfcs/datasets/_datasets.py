from pathlib import Path

HERE = Path(__file__).parent


def example() -> Path:
    """Example fcs file from FlowIO.

    https://github.com/whitews/FlowIO/blob/master/examples/fcs_files/100715.fcs
    """
    return HERE / "example.fcs"
