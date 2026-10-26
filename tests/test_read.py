import readfcs


def test_read_returns_writeable_X():
    adata = readfcs.read(readfcs.datasets.Oetjen18_t1())
    assert adata.X.flags.writeable
    assert adata.shape == (241552, 20)
    assert adata.obs_names[:3].tolist() == ["0", "1", "2"]
    assert adata.var_names[:3].tolist() == ["FSC-A", "FSC-H", "FSC-W"]
    # in-place operations on X must not fail
    adata.X[0, 0] = 1.0
