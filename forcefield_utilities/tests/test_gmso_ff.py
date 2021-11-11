import pytest

from forcefield_utilities.prepackaged import FoyerFFs
from forcefield_utilities.tests.base_test import BaseTest


class TestGMSOFFConversion(BaseTest):
    @pytest.fixture(scope="session")
    def oplsaa_gmso(self):
        return FoyerFFs.oplsaa.to_gmso_ff()

    @pytest.fixture(scope="session")
    def gaff_gmso(self):
        return FoyerFFs.gaff.to_gmso_ff()

    def test_atom_types(self):
        pass

    def test_bond_types(self):
        pass

    def test_angle_types(self):
        pass

    def test_dihedral_types(self):
        pass

    def test_improper_types(self):
        pass

    def test_metadata(self):
        pass
