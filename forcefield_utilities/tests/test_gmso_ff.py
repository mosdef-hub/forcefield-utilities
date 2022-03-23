import numpy as np
import pytest
import unyt as u
from gmso.utils._constants import FF_TOKENS_SEPARATOR

from forcefield_utilities.foyer_xml import AtomTypes, NonBondedForce
from forcefield_utilities.xml_loader import FoyerFFs
from forcefield_utilities.tests.base_test import BaseTest

parameters_map = {"length": "r_eq", "angle": "theta_eq"}


def assert_atomtypes_equivalency(parameters_ff, gmso_ff):
    atom_types_gmso = gmso_ff.atom_types
    non_bonded_forces = list(
        filter(lambda c: isinstance(c, NonBondedForce), parameters_ff.children)
    ).pop()
    xml_atom_types = list(
        filter(lambda c: isinstance(c, AtomTypes), parameters_ff.children)
    ).pop()
    nb_atom_types = {
        non_bonded.atom_type: non_bonded for non_bonded in non_bonded_forces.children
    }

    xml_atom_types = {
        xml_atom_type.name: xml_atom_type for xml_atom_type in xml_atom_types.children
    }
    for atom_type_name in nb_atom_types:
        type_ = nb_atom_types[atom_type_name]
        atom_type_gmso = atom_types_gmso[atom_type_name]
        assert type_.atom_type == atom_type_gmso.name
        assert np.allclose(
            type_.charge, (atom_type_gmso.charge / u.elementary_charge).value
        )
        assert type_.sigma == atom_type_gmso.parameters["sigma"].value
        assert type_.epsilon == atom_type_gmso.parameters["epsilon"].value

    for atom_type_name in xml_atom_types:
        atom_type = xml_atom_types[atom_type_name]
        atom_type_gmso = atom_types_gmso[atom_type_name]
        assert atom_type.name == atom_type_gmso.name
        assert atom_type.atom_class == atom_type_gmso.atomclass
        assert atom_type.doi == atom_type_gmso.doi
        assert atom_type.smarts_def == atom_type_gmso.definition
        assert atom_type.desc == atom_type_gmso.description
        assert atom_type.element == atom_type_gmso.get_tag("element")
        assert np.allclose(atom_type.mass, atom_type_gmso.mass.value)
        if atom_type.overrides:
            for splited in atom_type.overrides.split(","):
                assert splited.strip() in atom_type_gmso.overrides
        else:
            assert atom_type_gmso.overrides == set()


def assert_forces_equivalency(
    parameters_ff,
    gmso_ff,
    xml_force_attr="HarmonicBondForce",
    gmso_potential_type="bond_types",
    has_mixed_children=False,
):
    all_children_iter = parameters_ff.iterate_on(xml_force_attr)
    gmso_potentials = getattr(gmso_ff, gmso_potential_type)

    if len(gmso_potentials) == 0:
        return

    for child in all_children_iter:
        count = 0
        classes = []
        types = []
        while hasattr(child, f"class{count+1}"):
            classes.append(getattr(child, f"class{count+1}"))
            count += 1
        count = 0
        while hasattr(child, f"type{count + 1}"):
            types.append(getattr(child, f"type{count + 1}"))
            count += 1

        if any(classes):
            classes = FF_TOKENS_SEPARATOR.join(
                class_ if class_ else "*" for class_ in classes
            )
            try:
                gmso_potential = gmso_potentials[classes]
            except KeyError as e:
                if has_mixed_children:
                    continue
                else:
                    raise e
            assert gmso_potential.member_types is None
            assert gmso_potential.member_classes == tuple(
                classes.split(FF_TOKENS_SEPARATOR)
            )
        else:
            types = FF_TOKENS_SEPARATOR.join(type_ if type_ else "*" for type_ in types)
            try:
                gmso_potential = gmso_potentials[types]
            except KeyError as e:
                if has_mixed_children:
                    continue
                else:
                    raise e
            assert gmso_potential.member_classes is None
            assert gmso_potential.member_types == tuple(
                types.split(FF_TOKENS_SEPARATOR)
            )

        btype_params = gmso_potential.get_parameters()
        btype_params_xml = child.parameters()
        for key in btype_params_xml:
            if key in btype_params:
                assert np.allclose(btype_params[key].value, btype_params_xml[key])
            else:
                assert key in parameters_map
                new_key_gmso = parameters_map[key]
                assert np.allclose(
                    btype_params[new_key_gmso].value, btype_params_xml[key]
                )


class TestGMSOFFConversionOPLSAA(BaseTest):
    @pytest.fixture(scope="session")
    def oplsaa_gmso(self):
        return FoyerFFs.get_ff("oplsaa").to_gmso_ff()

    def test_atom_types(self, oplsaa_gmso):
        assert_atomtypes_equivalency(FoyerFFs.get_ff("oplsaa"), oplsaa_gmso)

    def test_bond_types(self, oplsaa_gmso):
        assert_forces_equivalency(FoyerFFs.get_ff("oplsaa"), oplsaa_gmso)

    def test_angle_types(self, oplsaa_gmso):
        assert_forces_equivalency(
            FoyerFFs.get_ff("oplsaa"), oplsaa_gmso, "HarmonicAngleForce", "angle_types"
        )

    def test_dihedral_types(self, oplsaa_gmso):
        assert_forces_equivalency(
            FoyerFFs.get_ff("oplsaa"), oplsaa_gmso, "RBTorsionForce", "dihedral_types"
        )

    def test_metadata(self, oplsaa_gmso):
        assert oplsaa_gmso.name == FoyerFFs.get_ff("oplsaa").name
        assert oplsaa_gmso.version == FoyerFFs.get_ff("oplsaa").version
        non_bonded_forces = list(
            filter(
                lambda c: isinstance(c, NonBondedForce),
                FoyerFFs.get_ff("oplsaa").children,
            )
        ).pop()
        scaling_factors = {
            "electrostatics14Scale": non_bonded_forces.coulomb14scale,
            "nonBonded14Scale": non_bonded_forces.lj14scale,
        }
        assert scaling_factors == oplsaa_gmso.scaling_factors


class TestGMSOFFConversionTRAPPEUA(BaseTest):
    @pytest.fixture(scope="session")
    def trappe_ua_gmso(self):
        return FoyerFFs.get_ff("trappe_ua").to_gmso_ff()

    def test_atom_types(self, trappe_ua_gmso):
        assert_atomtypes_equivalency(FoyerFFs.get_ff("trappe_ua"), trappe_ua_gmso)

    def test_bond_types(self, trappe_ua_gmso):
        assert_forces_equivalency(FoyerFFs.get_ff("trappe_ua"), trappe_ua_gmso)

    def test_angle_types(self, trappe_ua_gmso):
        assert_forces_equivalency(
            FoyerFFs.get_ff("trappe_ua"),
            trappe_ua_gmso,
            "HarmonicAngleForce",
            "angle_types",
        )

    def test_dihedral_types(self, trappe_ua_gmso):
        assert_forces_equivalency(
            FoyerFFs.get_ff("trappe_ua"),
            trappe_ua_gmso,
            "RBTorsionForce",
            "dihedral_types",
        )

    def test_metadata(self, trappe_ua_gmso):
        assert trappe_ua_gmso.name == FoyerFFs.get_ff("trappe_ua").name
        assert trappe_ua_gmso.version == FoyerFFs.get_ff("trappe_ua").version
        non_bonded_forces = list(
            filter(
                lambda c: isinstance(c, NonBondedForce),
                FoyerFFs.get_ff("trappe_ua").children,
            )
        ).pop()
        scaling_factors = {
            "electrostatics14Scale": non_bonded_forces.coulomb14scale,
            "nonBonded14Scale": non_bonded_forces.lj14scale,
        }
        assert scaling_factors == trappe_ua_gmso.scaling_factors


@pytest.mark.skipif(
    condition=not hasattr(FoyerFFs, "gaff"),
    reason="Gaff Forcefield not found",
)
class TestGMSOFFConversionGAFF:
    @pytest.fixture(scope="session")
    def gaff_gmso(self):
        return FoyerFFs.gaff.to_gmso_ff()

    def test_atom_types(self, gaff_gmso):
        assert_atomtypes_equivalency(FoyerFFs.gaff, gaff_gmso)

    def test_bond_types(self, gaff_gmso):
        assert_forces_equivalency(FoyerFFs.gaff, gaff_gmso)

    def test_angle_types(self, gaff_gmso):
        assert_forces_equivalency(
            FoyerFFs.gaff, gaff_gmso, "HarmonicAngleForce", "angle_types"
        )

    def test_dihedral_types(self, gaff_gmso):
        assert_forces_equivalency(
            FoyerFFs.get_ff('trappe_ua'),
            gaff_gmso,
            "PeriodicTorsionForce",
            "dihedral_types",
            has_mixed_children=True,
        )

    def test_dihedral_types(self, gaff_gmso):
        assert_forces_equivalency(
            FoyerFFs.get_ff('trappe_ua'),
            gaff_gmso,
            "PeriodicTorsionForce",
            "improper_types",
            has_mixed_children=True,
        )

    def test_metadata(self, gaff_gmso):
        assert gaff_gmso.name == FoyerFFs.gaff.name
        assert gaff_gmso.version == FoyerFFs.gaff.version
        non_bonded_forces = list(
            filter(lambda c: isinstance(c, NonBondedForce), FoyerFFs.gaff.children)
        ).pop()
        scaling_factors = {
            "electrostatics14Scale": non_bonded_forces.coulomb14scale,
            "nonBonded14Scale": non_bonded_forces.lj14scale,
        }
        assert scaling_factors == gaff_gmso.scaling_factors
