import pytest
import unyt as u
from gmso.tests.utils import get_path
from lxml import etree

from forcefield_utilities.gmso_xml import ForceField
from forcefield_utilities.tests.base_test import BaseTest


class TestGMSOXML(BaseTest):
    @pytest.fixture(scope="session")
    def ff_example_zero(self):
        example_zero = get_path("ethylene.xml")
        with open(example_zero, "r") as example_zero_file:
            root = etree.parse(example_zero_file).getroot()
        return ForceField.load_from_etree(root).to_gmso_ff()

    def test_metadata(self, ff_example_zero):
        assert ff_example_zero.scaling_factors == {
            "electrostatic_14_scale": 0.5,
            "nonbonded_14_scale": 0.5,
        }

    def test_atom_types(self, ff_example_zero):
        opls_143 = ff_example_zero.atom_types["opls_143"]
        assert opls_143.name == "opls_143"
        assert opls_143.charge == -0.23 * u.elementary_charge
        assert opls_143.get_tag("element") == "C"
