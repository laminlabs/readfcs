import flowio
import pytest

import readfcs


def test_read_returns_writeable_X():
    adata = readfcs.read(readfcs.datasets.Oetjen18_t1())
    assert adata.X.flags.writeable
    assert adata.shape == (241552, 20)
    assert adata.obs_names[:3].tolist() == ["0", "1", "2"]
    assert adata.var_names[:3].tolist() == ["FSC-A", "FSC-H", "FSC-W"]
    # in-place operations on X must not fail
    adata.X[0, 0] = 1.0


def test_read_anndata_writes_h5ad(tmp_path):
    adata = readfcs.read(readfcs.datasets.Oetjen18_t1())
    spill = adata.uns["meta"]["spill"]
    assert "LIVE_DEAD" in spill.columns
    assert not any("/" in str(label) for label in spill.columns)
    assert not any("/" in str(label) for label in spill.index)
    path = tmp_path / "oetjen18_t1.h5ad"
    adata.write(path)
    assert path.exists()



@pytest.mark.parametrize(
    ("reindex", "expected_index"),
    [(True, ["FSC-A", "SSC-A"]), (False, ["1", "2"])],
)
def test_read_without_marker_names(tmp_path, reindex, expected_index):
    path = tmp_path / "no_markers.fcs"
    channels = ["FSC-A", "SSC-A"]
    with path.open("wb") as file:
        flowio.create_fcs(file, [1.0, 2.0, 3.0, 4.0], channels)

    adata = readfcs.read(path, reindex=reindex)

    assert adata.X.tolist() == [[1.0, 2.0], [3.0, 4.0]]
    assert adata.var["channel"].tolist() == channels
    assert adata.var["marker"].tolist() == ["", ""]
    assert adata.var_names.tolist() == expected_index
