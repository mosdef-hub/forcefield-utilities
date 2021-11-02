import foyer
import pytest


class BaseTest:
    @pytest.fixture(scope="session")
    def oplsaa_foyer(self):
        return foyer.forcefields.load_OPLSAA()

    @pytest.fixture(scope="session")
    def gaff_foyer(self):
        return foyer.forcefields.load_GAFF()

    @pytest.fixture(scope="session")
    def missing_atom_classes_oplsaa(self):
        return {
            "NZ",
            "H",
            "C3",
            "CM",
            "CR",
            "H3",
            "OW",
            "OU",
            "CT_3",
            "HA",
            "CT_4",
            "OY",
            "CK",
            "CU",
            "CA",
            "C+",
            "O2",
            "CV",
            "N2",
            "LP",
            "CT_2",
            "CW",
            "CQ",
            "N",
            "C*",
            "Br",
            "C",
            "H4",
            "HC",
            "NY",
            "CB",
            "CT",
            "O",
            "C_3",
            "SY",
            "U",
            "C=",
            "N3",
            "P",
            "CN",
            "OH",
            "S",
            "NB",
            "CX",
            "N*",
        }
