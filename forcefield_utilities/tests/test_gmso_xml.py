import pytest
import unyt as u
from gmso.tests.utils import get_path
from lxml import etree
from sympy import sympify

from forcefield_utilities.tests.base_test import BaseTest
from forcefield_utilities.tests.utils import get_test_file_path
from forcefield_utilities.xml_loader import GMSOFFs


class TestEthyleneFF(BaseTest):
    @pytest.fixture(scope="session")
    def ff_example_zero(self):
        example_zero = get_path("ethylene.xml")
        return GMSOFFs().load(example_zero).to_gmso_ff()

    def test_metadata(self, ff_example_zero):
        assert ff_example_zero.scaling_factors == {
            "electrostatics14Scale": 0.5,
            "nonBonded14Scale": 0.5,
        }
        assert ff_example_zero.combining_rule == "geometric"

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


class TestTwoPropanolMIEFF(BaseTest):
    @pytest.fixture(scope="session")
    def propanol_ua_mie(self):
        propanol_ua_mie_path = get_test_file_path("propanol_Mie_ua.xml")
        return GMSOFFs().load(propanol_ua_mie_path).to_gmso_ff()

    def test_metadata(self, propanol_ua_mie):
        assert (
            propanol_ua_mie.name
            == "Mie two-propanol- This is for testing only and not for use for simulations "
        )
        assert propanol_ua_mie.scaling_factors == {
            "electrostatics14Scale": 0.0,
            "nonBonded14Scale": 0.0,
        }
        assert propanol_ua_mie.combining_rule == "geometric"

    def test_atom_types(self, propanol_ua_mie):
        ch3_sp3 = propanol_ua_mie.atom_types["CH3_sp3"]

        assert ch3_sp3.name == "CH3_sp3"
        assert ch3_sp3.atomclass == "CH3"
        assert ch3_sp3.get_tag("element") == "_CH3"
        assert u.allclose_units(ch3_sp3.charge, 0.0 * u.coulomb)
        assert ch3_sp3.definition == "[_CH3;X1][_CH3,_HC]"
        assert u.allclose_units(ch3_sp3.mass, 15.03500 * u.amu)
        assert (
            ch3_sp3.description
            == "Alkane CH3, Mie using the k constant from Trappe-UA"
        )
        assert ch3_sp3.doi == "10.1021/jp984742e and 10.1021/jp972543+"
        assert ch3_sp3.overrides == set()

        assert ch3_sp3.expression == sympify(
            "(n/(n-m)) * (n/m)**(m/(n-m)) * epsilon * ((sigma/r)**n - (sigma/r)**m)"
        )

        parameters = ch3_sp3.get_parameters()
        assert u.allclose_units(
            parameters["epsilon"], 0.194746017346801 * u.kcal / u.mol
        )
        assert u.allclose_units(parameters["sigma"], 3.751 * u.angstrom)
        assert u.allclose_units(parameters["n"], 11 * u.dimensionless)
        assert u.allclose_units(parameters["m"], 6 * u.dimensionless)

        o = propanol_ua_mie.atom_types["O"]

        assert o.name == "O"
        assert o.atomclass == "O"
        assert o.get_tag("element") == "O"
        assert u.allclose_units(o.charge, -0.700 * u.elementary_charge)
        assert o.definition == "OH"
        assert u.allclose_units(o.mass, 15.99940 * u.amu)
        assert o.description == "Oxygen in hydroxyl"
        assert o.doi == "10.1021/jp003882x"
        assert o.overrides == set()

        assert o.expression == sympify(
            "(n/(n-m)) * (n/m)**(m/(n-m)) * epsilon * ((sigma/r)**n - (sigma/r)**m)"
        )

        parameters = o.get_parameters()
        assert u.allclose_units(
            parameters["epsilon"], 0.184809996053596 * u.kcal / u.mol
        )
        assert u.allclose_units(parameters["sigma"], 3.021 * u.angstrom)
        assert u.allclose_units(parameters["n"], 13 * u.dimensionless)
        assert u.allclose_units(parameters["m"], 6 * u.dimensionless)

    def test_bond_types(self, propanol_ua_mie):
        bond_type_ch3_ch = propanol_ua_mie.bond_types["CH3~CH"]
        assert bond_type_ch3_ch.name == "BondType_Harmonic_CH3_CH"
        assert bond_type_ch3_ch.member_classes == ("CH3", "CH")
        assert bond_type_ch3_ch.expression == sympify("k * (r-r_eq)**2")

        parameters = bond_type_ch3_ch.get_parameters()
        assert u.allclose_units(
            parameters["k"], 1200.80305927342 * u.kcal / u.mol / u.angstrom**2
        )
        assert u.allclose_units(parameters["r_eq"], 1.5401 * u.angstrom)

    def test_angle_types(self, propanol_ua_mie):
        angle_type_ch3_ch_o = propanol_ua_mie.angle_types["CH3~CH~O"]

        assert angle_type_ch3_ch_o.name == "AngleType_Harmonic_CH3_CH_O"
        assert angle_type_ch3_ch_o.member_classes == ("CH3", "CH", "O")
        assert angle_type_ch3_ch_o.expression == sympify(
            "k * (theta - theta_eq)**2"
        )

        parameters = angle_type_ch3_ch_o.get_parameters()
        assert u.allclose_units(
            parameters["k"], 100.155094635497 * u.kcal / u.mol / u.radian**2
        )
        assert u.allclose_units(parameters["theta_eq"], 109.51 * u.degree)

    def test_dihedral_types(self, propanol_ua_mie):
        dihedral_type_ch3_ch_o_h = propanol_ua_mie.dihedral_types["CH3~CH~O~H"]

        assert (
            dihedral_type_ch3_ch_o_h.name
            == "DihedralType_Periodic_Proper_CH3_CH_O_H"
        )
        assert dihedral_type_ch3_ch_o_h.member_classes == (
            "CH3",
            "CH",
            "O",
            "H",
        )
        assert dihedral_type_ch3_ch_o_h.expression == sympify(
            "k0 + k1 * (1 + cos(1 * phi - phi_eq1)) + "
            "k2 * (1 + cos(2 * phi - phi_eq2)) + "
            "k3 * (1 + cos(3 * phi - phi_eq3)) + "
            "k4 * (1 + cos(4 * phi - phi_eq4)) + "
            "k5 * (1 + cos(5 * phi - phi_eq5))"
        )

        parameters = dihedral_type_ch3_ch_o_h.get_parameters()
        assert u.allclose_units(parameters["k0"], 0.0 * u.kcal / u.mol)
        assert u.allclose_units(
            parameters["k1"], 0.416955197548017 * u.kcal / u.mol
        )
        assert u.allclose_units(
            parameters["k2"], -0.0579667482245528 * u.kcal / u.mol
        )
        assert u.allclose_units(
            parameters["k3"], 0.37345529632637 * u.kcal / u.mol
        )
        assert u.allclose_units(parameters["k4"], 0.0 * u.kcal / u.mol)
        assert u.allclose_units(parameters["k5"], 0.0 * u.kcal / u.mol)

        assert u.allclose_units(parameters["phi_eq1"], 0.0 * u.degree)
        assert u.allclose_units(parameters["phi_eq2"], 180 * u.degree)
        assert u.allclose_units(parameters["phi_eq3"], 0.0 * u.degree)
        assert u.allclose_units(parameters["phi_eq4"], 0.0 * u.degree)
        assert u.allclose_units(parameters["phi_eq5"], 0.0 * u.degree)


class TestListParameters(BaseTest):
    @pytest.fixture(scope="session")
    def propanol_ua_mie_list(self):
        propanol_ua_mie_path = get_test_file_path(
            "propanol_Mie_ua_list_params.xml"
        )
        return GMSOFFs().load(propanol_ua_mie_path).to_gmso_ff()

    def test_dihedral_params(self, propanol_ua_mie_list):
        dih_with_list = propanol_ua_mie_list.dihedral_types["CH3~CH~O~H"]
        params = dih_with_list.get_parameters()
        print(params)
        assert u.allclose_units(
            params["phi_eq"], [0.0, 180.0, 0.0, 0.0, 0.0, 0.0] * u.degree
        )
        assert u.allclose_units(
            params["k"],
            [0.0, 0.4169552, -0.05796675, 0.3734553, 0.0, 0.0] * u.kCal / u.mol,
        )
