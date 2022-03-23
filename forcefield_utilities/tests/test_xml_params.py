import pytest

from forcefield_utilities import FoyerFFs
from forcefield_utilities.foyer_xml import PeriodicProper, RBProper
from forcefield_utilities.tests.base_test import BaseTest


class TestXMLParams(BaseTest):
    def test_atom_types_oplsaa(self, oplsaa_foyer):
        for atom_type in FoyerFFs.get_ff("oplsaa").iterate_on(
            children_type="NonbondedForce"
        ):
            assert atom_type.dict(
                exclude_none=True, exclude={"atom_type"}
            ) == oplsaa_foyer.get_parameters(group="atoms", key=atom_type.atom_type)

    def test_bond_types_oplsaa(self, oplsaa_foyer, missing_atom_classes_oplsaa):
        for harmonic_bond in FoyerFFs.get_ff("oplsaa").iterate_on(
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
        for harmonic_angle in FoyerFFs.get_ff("oplsaa").iterate_on(
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

    def test_rb_torsions(self, oplsaa_foyer, missing_atom_classes_oplsaa):
        for rb_torsion in FoyerFFs.get_ff("oplsaa").iterate_on(
            children_type="RBTorsionForce"
        ):
            if (
                rb_torsion.class1 not in missing_atom_classes_oplsaa
                and rb_torsion.class2 not in missing_atom_classes_oplsaa
                and rb_torsion.class3 not in missing_atom_classes_oplsaa
                and rb_torsion.class4 not in missing_atom_classes_oplsaa
            ):
                ff_params = oplsaa_foyer.get_parameters(
                    group="rb_propers"
                    if isinstance(rb_torsion, RBProper)
                    else "rb_impropers",
                    key=[
                        rb_torsion.class1,
                        rb_torsion.class2,
                        rb_torsion.class3,
                        rb_torsion.class4,
                    ],
                    keys_are_atom_classes=True,
                )
                xml_params = rb_torsion.dict(
                    include={"c0", "c1", "c2", "c3", "c4", "c5"}
                )
                assert xml_params == ff_params

    @pytest.mark.skipif(
        condition=not hasattr(FoyerFFs, "gaff"),
        reason="Gaff Forcefield not found",
    )
    def test_atomtypes_gaff(self, gaff_foyer):
        for atom_type in FoyerFFs.gaff.iterate_on(children_type="NonbondedForce"):
            assert atom_type.dict(
                exclude_none=True, exclude={"atom_type"}
            ) == gaff_foyer.get_parameters(group="atoms", key=atom_type.atom_type)

    @pytest.mark.skipif(
        condition=not hasattr(FoyerFFs, "gaff"),
        reason="Gaff Forcefield not found",
    )
    def test_bond_types_gaff(self, gaff_foyer):
        for harmonic_bond in FoyerFFs.gaff.iterate_on(
            children_type="HarmonicBondForce"
        ):

            assert harmonic_bond.dict(
                exclude={"class1", "class2"}, exclude_none=True
            ) == gaff_foyer.get_parameters(
                group="harmonic_bonds",
                key=[harmonic_bond.class1, harmonic_bond.class2],
                keys_are_atom_classes=True,
            )

    @pytest.mark.skipif(
        condition=not hasattr(FoyerFFs, "gaff"),
        reason="Gaff Forcefield not found",
    )
    def test_angle_types_gaff(self, gaff_foyer):
        for harmonic_angle in FoyerFFs.gaff.iterate_on(
            children_type="HarmonicAngleForce"
        ):
            angle_params = harmonic_angle.dict(
                exclude={"class1", "class2", "class3"}, exclude_none=True
            )
            angle_params["theta"] = angle_params.pop("angle")

            assert angle_params == gaff_foyer.get_parameters(
                group="harmonic_angles",
                key=[
                    harmonic_angle.class1,
                    harmonic_angle.class2,
                    harmonic_angle.class3,
                ],
                keys_are_atom_classes=True,
            )

    # FixME: Unconditionally Skipped
    # @pytest.mark.skipif(
    #     condition=not hasattr(FoyerFFs, "gaff"),
    #     reason="Gaff Forcefield not found",
    # )
    @pytest.mark.skip
    def test_periodic_torsions_gaff(self, gaff_foyer):
        for periodic_torsion in FoyerFFs.gaff.iterate_on(
            children_type="PeriodicTorsionForce"
        ):
            key = [
                periodic_torsion.class1,
                periodic_torsion.class2,
                periodic_torsion.class3,
                periodic_torsion.class4,
            ]
            if key != ["", "c", "c1", ""]:
                ff_params = gaff_foyer.get_parameters(
                    group="periodic_propers"
                    if isinstance(periodic_torsion, PeriodicProper)
                    else "periodic_impropers",
                    key=key,
                    keys_are_atom_classes=True,
                )
                xml_params = periodic_torsion.dict(
                    include={"periodicity", "phase", "k"}
                )
                assert xml_params == ff_params
