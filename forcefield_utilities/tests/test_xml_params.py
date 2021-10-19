import foyer.exceptions

from forcefield_utilities import FoyerFFs
from forcefield_utilities.tests.base_test import BaseTest


class TestXMLParams(BaseTest):
    def test_atom_types_oplsaa(self, oplsaa_foyer):
        for atom_type in FoyerFFs.oplsaa.iterate_on(
            children_type="NonbondedForce"
        ):
            assert atom_type.dict(
                exclude_none=True, exclude={"atom_type"}
            ) == oplsaa_foyer.get_parameters(
                group="atoms", key=atom_type.atom_type
            )

    def test_bond_types_oplsaa(self, oplsaa_foyer, missing_atom_classes_oplsaa):
        for harmonic_bond in FoyerFFs.oplsaa.iterate_on(
            children_type="HarmonicBondForce"
        ):
            if (
                harmonic_bond.class1 not in missing_atom_classes_oplsaa
                and harmonic_bond.class2 not in missing_atom_classes_oplsaa
            ):
                assert harmonic_bond.dict(
                    exclude={"class1", "class2"}, exclude_none=True
                ) == oplsaa_foyer.get_parameters(
                    group="harmonic_bonds",
                    key=[harmonic_bond.class1, harmonic_bond.class2],
                    keys_are_atom_classes=True,
                )

    def test_harmonic_angle_types_oplsaa(
        self, oplsaa_foyer, missing_atom_classes_oplsaa
    ):
        for harmonic_angle in FoyerFFs.oplsaa.iterate_on(
            children_type="HarmonicAngleForce"
        ):
            if (
                harmonic_angle.class1 not in missing_atom_classes_oplsaa
                and harmonic_angle.class2 not in missing_atom_classes_oplsaa
                and harmonic_angle.class3 not in missing_atom_classes_oplsaa
            ):
                angle_params = harmonic_angle.dict(
                    exclude={"class1", "class2", "class3"}, exclude_none=True
                )
                angle_params["theta"] = angle_params.pop("angle")
                assert angle_params == oplsaa_foyer.get_parameters(
                    group="harmonic_angles",
                    key=[
                        harmonic_angle.class1,
                        harmonic_angle.class2,
                        harmonic_angle.class3,
                    ],
                    keys_are_atom_classes=True,
                )
