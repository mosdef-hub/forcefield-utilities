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
        angle_type_harmonic_2 = ff_example_zero.angle_types[
            "opls_144~opls_143~opls_144"
        ]
        assert (
            str(angle_type_harmonic_2.expression)
            == "0.5*k*(theta - theta_eq)**2"
        )
        assert angle_type_harmonic_2.member_types == (
            "opls_144",
            "opls_143",
            "opls_144",
        )
        assert angle_type_harmonic_2.parameters == {
            "k": u.unyt_quantity(292.88, "kJ/(mol*radian**2)"),
            "theta_eq": u.unyt_quantity(2.0420352248, "radian"),
        }

        assert angle_type_harmonic_2.name == "AngleType-Harmonic-2"

    def test_dihedral_types(self, ff_example_zero):
        dihedral_type_rb_1 = ff_example_zero.dihedral_types[
            "opls_144~opls_143~opls_143~opls_144"
        ]

        assert (
            str(dihedral_type_rb_1.expression)
            == "c_0 + c_1*cos(psi) + c_2*cos(psi)**2 + c_3*cos(psi)**3 + c_4*cos(psi)**4 + c_5*cos(psi)**5"
        )
        assert dihedral_type_rb_1.member_types == (
            "opls_144",
            "opls_143",
            "opls_143",
            "opls_144",
        )

        assert dihedral_type_rb_1.parameters == {
            "c_0": u.unyt_quantity(58.576, "kJ/mol"),
            "c_1": u.unyt_quantity(0.0, "kJ/mol"),
            "c_2": u.unyt_quantity(-58.576, "kJ/mol"),
            "c_3": u.unyt_quantity(0.0, "kJ/mol"),
            "c_4": u.unyt_quantity(0.0, "kJ/mol"),
            "c_5": u.unyt_quantity(0.0, "kJ/mol"),
        }

        assert dihedral_type_rb_1.name == "DihedralType-RyckaertBellemans-1"

    def test_improper_types(self, ff_example_zero):
        pass
