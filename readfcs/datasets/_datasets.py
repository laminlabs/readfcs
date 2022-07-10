from pathlib import Path

HERE = Path(__file__).parent


def example() -> Path:
    return HERE / "example.fcs"
