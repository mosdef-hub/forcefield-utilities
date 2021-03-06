import pytest

from forcefield_utilities.tests.base_test import BaseTest
from forcefield_utilities.utils import get_package_file_path
from forcefield_utilities.xml_loader import FoyerFFs


class TestXMLLoader(BaseTest):
    @pytest.fixture
    def foyer_xml_loader(self):
        return FoyerFFs()

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
