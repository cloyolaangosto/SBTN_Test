import pytest

from sbtn_leaf import data_loader


@pytest.fixture(autouse=True)
def clear_loader_caches():
    """Ensure data loader caches are reset between tests."""

    data_loader._load_crop_naming_index_table.cache_clear()
    data_loader._load_fao_crop_yields_table.cache_clear()
    data_loader._load_ecoregions_shapefile.cache_clear()

    yield

    data_loader._load_crop_naming_index_table.cache_clear()
    data_loader._load_fao_crop_yields_table.cache_clear()
    data_loader._load_ecoregions_shapefile.cache_clear()


def test_crop_naming_index_missing_file(monkeypatch, tmp_path):
    monkeypatch.setattr(data_loader, "DATA_DIR", tmp_path)
    with pytest.raises(FileNotFoundError) as exc:
        data_loader.get_crop_naming_index_table()
    message = str(exc.value)
    assert "crop naming index table" in message
    assert str(tmp_path) in message


def test_ecoregions_shapefile_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(data_loader, "DATA_DIR", tmp_path)
    with pytest.raises(FileNotFoundError) as exc:
        data_loader.get_ecoregions_shapefile()
    message = str(exc.value)
    assert "ecoregions shapefile" in message
    assert str(tmp_path) in message
