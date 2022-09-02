import pytest

from foyer.forcefield import get_available_forcefield_loaders

from forcefield_utilities.tests.base_test import BaseTest
from forcefield_utilities.tests.utils import get_test_file_path
from forcefield_utilities.utils import get_package_file_path
from forcefield_utilities.xml_loader import FoyerFFs, GMSOFFs


class TestXMLLoader(BaseTest):
    @pytest.fixture
    def foyer_xml_loader(self):
        return FoyerFFs()

    @pytest.fixture
    def gmso_xml_loader(self):
        return GMSOFFs()

    @pytest.mark.skipif(
    condition="load_GAFF"
    not in map(lambda func: func.__name__, get_available_forcefield_loaders()),
    reason="GAFF Plugin is not installed",
    )
    def test_load_gaff(self, foyer_xml_loader):
        foyer_xml_loader.get_ff("gaff")
        assert "gaff" in foyer_xml_loader.loaded_ffs

    def test_caching(self, foyer_xml_loader):
        oplsaa = foyer_xml_loader.load("oplsaa")
        id_oplsaa = id(oplsaa)
        for j in range(100):
            assert id_oplsaa == id(foyer_xml_loader.load("oplsaa"))
        foyer_xml_loader.clear_cache()
        assert foyer_xml_loader.loaded_ffs == {}
        assert id_oplsaa != foyer_xml_loader.load("oplsaa")

    def test_getitem_get_attr(self, foyer_xml_loader):
        oplsaa = foyer_xml_loader.load("oplsaa")
        assert (
            id(oplsaa)
            == id(foyer_xml_loader.oplsaa)
            == id(foyer_xml_loader["oplsaa"])
        )
        foyer_xml_loader.clear_cache()
        with pytest.raises(AttributeError):
            foyer_xml_loader.oplsaa
        with pytest.raises(KeyError):
            foyer_xml_loader["oplsaa"]

    def test_deprecation_warning(self, foyer_xml_loader):
        with pytest.warns(DeprecationWarning):
            foyer_xml_loader.load("trappe-ua", rel_to_module=True)

    def test_custom_register(self, foyer_xml_loader):
        xml_path = get_package_file_path("foyer", "tests/files/benzene_lb.xml")
        foyer_xml_loader.register_custom_forcefield("benzene_lb", xml_path)
        foyer_xml_loader.load("benzene_lb")
        assert "benzene_lb" in foyer_xml_loader.loaded_ffs
        with pytest.raises(ValueError):
            foyer_xml_loader.register_custom_forcefield(
                "benzene_lb", xml_path, overwrite=False
            )
        loaded_id = id(foyer_xml_loader.load("benzene_lb"))
        foyer_xml_loader.register_custom_forcefield(
            "benzene_lb", xml_path, overwrite=True
        )

        assert "benzene_lb" in foyer_xml_loader.loaded_ffs
        assert loaded_id != id(foyer_xml_loader.load("benzene_lb"))

    def test_different_loading_entries(self, foyer_xml_loader, gmso_xml_loader):
        assert id(foyer_xml_loader.loaded_ffs) != id(gmso_xml_loader.loaded_ffs)
        gmso_xml_loader.load(get_test_file_path("propanol_Mie_ua.xml"))
        assert "propanol_Mie_ua" in gmso_xml_loader.loaded_ffs
        assert "propanol_Mie_ua" not in foyer_xml_loader.loaded_ffs

    def test_class_methods_gmso_ff(self):
        gmso_xml_loader = GMSOFFs()
        gmso_xml_loader.load(get_test_file_path("propanol_Mie_ua.xml"))
        ff1 = gmso_xml_loader.get_ff("propanol_Mie_ua")
        ff2 = gmso_xml_loader.load("propanol_Mie_ua")
        assert id(ff1) == id(ff2)
