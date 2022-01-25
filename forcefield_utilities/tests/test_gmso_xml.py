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

    def test_bond_types(self, ff_example_zero):
        btype_harmonic_1 = ff_example_zero.bond_types["opls_143~opls_143"]
        assert str(btype_harmonic_1.expression) == "0.5*k*(r - r_eq)**2"
        assert btype_harmonic_1.member_types == ("opls_143", "opls_143")
        assert btype_harmonic_1.parameters == {
            "k": u.unyt_quantity(459403.2, "kJ/mol/nm**2"),
            "r_eq": u.unyt_quantity(0.134, "nm"),
        }

        assert btype_harmonic_1.name == "BondType-Harmonic-1"

    def test_angle_types(self, ff_example_zero):
        pass
