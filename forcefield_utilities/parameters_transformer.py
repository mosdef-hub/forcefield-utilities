import unyt as u


class TransformationError(Exception):
    pass


class ParametersTransformer:
    @staticmethod
    def transform(name, parameters):
        if name == "LennardJonesPotential":
            return ParametersTransformer.transform_lennard_jones(parameters)
        elif name == "HarmonicBondPotential":
            return ParametersTransformer.transform_harmonic_bond(parameters)
        elif name == "HarmonicAnglePotential":
            return ParametersTransformer.transform_harmonic_angle(parameters)
        elif name == "RyckaertBellemansTorsionPotential":
            return ParametersTransformer.transform_rb_torsion(parameters)
        elif name == "PeriodicTorsionPotential":
            return ParametersTransformer.transform_periodic_torsion(parameters)
        else:
            raise TransformationError(
                f"No transformation is defined for {name}"
            )

    @staticmethod
    def transform_lennard_jones(parameters):
        transformed = {
            "epsilon": parameters["epsilon"] * u.kJ / u.mol,
            "sigma": parameters["sigma"] * u.nm,
        }
        return transformed

    @staticmethod
    def transform_harmonic_bond(parameters):
        transformed = {
            "k": parameters["k"] * u.kJ / (u.nm**2) / u.mol,
            "r_eq": parameters["length"] * u.nm,
        }
        return transformed

    @staticmethod
    def transform_harmonic_angle(parameters):
        transformed = {
            "k": parameters["k"] * u.kJ / (u.radian**2) / u.mol,
            "theta_eq": parameters["angle"] * u.radian,
        }
        return transformed

    @staticmethod
    def transform_rb_torsion(parameters):
        transformed = {
            "c0": parameters["c0"] * u.kJ / u.mol,
            "c1": parameters["c1"] * u.kJ / u.mol,
            "c2": parameters["c2"] * u.kJ / u.mol,
            "c3": parameters["c3"] * u.kJ / u.mol,
            "c4": parameters["c4"] * u.kJ / u.mol,
            "c5": parameters["c5"] * u.kJ / u.mol,
        }
        return transformed

    @staticmethod
    def transform_periodic_torsion(parameters):
        transformed = {
            "k": parameters["k"] * u.kJ / u.mol,
            "phi_eq": parameters["phase"] * u.radian,
            "n": parameters["periodicity"] * u.dimensionless,
        }
        return transformed
