from functools import lru_cache
from typing import List, Optional, Set, Union

import numpy as np
import sympy
import unyt as u
from gmso import ForceField as GMSOForceField
from gmso.core.angle_type import AngleType as GMSOAngleType
from gmso.core.atom_type import AtomType as GMSOAtomType
from gmso.core.bond_type import BondType as GMSOBondType
from gmso.core.dihedral_type import DihedralType as GMSODihedralType
from gmso.core.improper_type import ImproperType as GMSOImproperType
from gmso.utils._constants import FF_TOKENS_SEPARATOR
from pydantic import BaseModel, Field

from forcefield_utilities.utils import pad_with_wildcards


@lru_cache(maxsize=128)
def indep_vars(expr: str, dependent: frozenset) -> Set:
    """Given an expression and dependent variables, return independent variables for it."""
    return sympy.sympify(expr).free_symbols - dependent


class GMSOXMLTag(BaseModel):
    def parameters(self, units=None):
        params = self.children[0]
        params_dict = {}
        for parameter in params.children:
            if units is None:
                params_dict[parameter.name] = parameter.value
            else:
                params_dict[parameter.name] = (
                    parameter.value * units[parameter.name]
                )
        return params_dict

    class Config:
        arbitrary_types_allowed = True
        allow_population_by_field_name = True


class GMSOXMLChild(GMSOXMLTag):
    pass


class ParametersUnitDef(GMSOXMLTag):
    parameter: str = Field(
        ..., description="The name of the parameter", alias="parameter"
    )

    unit: str = Field(
        ..., description="The unit of the parameter", alias="unit"
    )

    @classmethod
    def load_from_etree(cls, root):
        return cls(**root.attrib)


class Parameter(GMSOXMLTag):
    name: str = Field(
        ..., description="The name of the parameter", alias="name"
    )

    value: Union[float, np.ndarray] = Field(
        ..., description="The value of the parameter", alias="value"
    )

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        if "value" in root.attrib:
            return cls(**attribs)
        else:
            children = root.getchildren()
            if len(children) == 0:
                raise ValueError(
                    "Neither a single value nor a sequence of values provided for "
                    f"parameter {attribs['name']}. Please provide one or the other"
                )
            value = np.array(
                [param_value.text for param_value in children], dtype=float
            )
            return cls(name=attribs["name"], value=value)


class Parameters(GMSOXMLTag):
    children: List[Parameter] = Field(
        ..., description="Parameter", alias="parameter"
    )

    @classmethod
    def load_from_etree(cls, root):
        children = []
        for el in root.iterchildren():
            if el.tag == "Parameter":
                children.append(Parameter.load_from_etree(el))
        return cls(children=children)


class AtomType(GMSOXMLTag):
    name: str = Field(
        ..., description="The name for this atom type", alias="name"
    )

    element: Optional[str] = Field(
        None, description="The element of the atom type", alias="element"
    )

    charge: Optional[float] = Field(
        None, description="The charge of the atom type", alias="charge"
    )

    mass: Optional[float] = Field(
        None, description="The mass of the atom type", alias="mass"
    )

    expression: Optional[str] = Field(
        None,
        description="The expression for this atom type",
        alias="expression",
    )

    independent_variables: Optional[str] = Field(
        None,
        description="The independent variables for this atom type",
        alias="independent_variables",
    )

    atomclass: Optional[str] = Field(
        None, description="The atomclass of this atomtype", alias="atomclass"
    )

    doi: Optional[str] = Field(
        None, description="The doi of this atomtype", alias="doi"
    )

    overrides: Optional[str] = Field(
        None, description="The overrides", alias="overrides"
    )

    definition: Optional[str] = Field(
        None,
        description="The smarts definition of this atom type",
        alias="definition",
    )

    description: Optional[str] = Field(
        None,
        description="The description of this atom type",
        alias="description",
    )

    children: List[Parameters] = Field(
        ..., description="The parameters and their values", alias="children"
    )

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "Parameters":
                children.append(Parameters.load_from_etree(el))
        return cls(children=children, **attribs)


class AtomTypes(GMSOXMLChild):
    name: Optional[str] = Field(
        None, description="The name for this atom type group", alias="name"
    )

    expression: str = Field(
        ...,
        description="The expression for this atom type group",
        alias="expression",
    )

    children: List[Union[ParametersUnitDef, AtomType]] = Field(
        ..., description="The children of AtomTypes", alias="parameters"
    )

    def to_gmso_potentials(self, default_units):
        potentials = {"atom_types": {}}
        parameters_units = filter(
            lambda c: isinstance(c, ParametersUnitDef), self.children
        )
        units = {
            parameter_unit.parameter: u.Unit(parameter_unit.unit)
            for parameter_unit in parameters_units
        }

        for atom_type in filter(
            lambda c: isinstance(c, AtomType), self.children
        ):
            atom_type_dict = atom_type.dict(
                by_alias=True,
                exclude={"children", "element"},
                exclude_none=True,
            )

            overrides = atom_type_dict.get("overrides")
            if overrides:
                atom_type_dict["overrides"] = set(
                    o.strip() for o in overrides.split(",")
                )
            else:
                atom_type_dict["overrides"] = set()

            if "expression" not in atom_type_dict:
                atom_type_dict["expression"] = self.expression
            atom_type_dict["parameters"] = atom_type.parameters(units)

            if not atom_type_dict.get("independent_variables"):
                atom_type_dict["independent_variables"] = indep_vars(
                    atom_type_dict["expression"],
                    frozenset(atom_type_dict["parameters"]),
                )

            if default_units.get("charge") and atom_type_dict.get("charge"):
                atom_type_dict["charge"] = (
                    atom_type_dict["charge"] * default_units["charge"]
                )

            if default_units.get("mass") and atom_type_dict.get("mass"):
                atom_type_dict["mass"] = (
                    atom_type_dict["mass"] * default_units["mass"]
                )
            gmso_atom_type = GMSOAtomType(**atom_type_dict)
            element = atom_type.element
            if element:
                gmso_atom_type.add_tag("element", element)
            potentials["atom_types"][atom_type.name] = gmso_atom_type
        return potentials

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "AtomType":
                children.append(AtomType.load_from_etree(el))
        return cls(children=children, **attribs)


class BondType(GMSOXMLTag):
    name: str = Field(
        None, description="The name of the bond type", alias="name"
    )

    class1: Optional[str] = Field(
        None, description="Class 1 for this bond type", alias="class1"
    )

    class2: Optional[str] = Field(
        None, description="Class 2 for this bond type", alias="class2"
    )

    type1: Optional[str] = Field(
        None, description="Type 1 for this bond type", alias="type1"
    )

    type2: Optional[str] = Field(
        None, description="Type 2 for this bond type", alias="type2"
    )

    children: List[Parameters] = Field(
        ..., description="The parameters and their values", alias="children"
    )

    @classmethod
    def load_from_etree(cls, root):
        children = []
        attribs = pad_with_wildcards(root.attrib, 2)
        for el in root.iterchildren():
            if el.tag == "Parameters":
                children.append(Parameters.load_from_etree(el))
        return cls(children=children, **attribs)


class BondTypes(GMSOXMLChild):
    name: Optional[str] = Field(
        None, description="The name for this bond types group", alias="name"
    )

    expression: str = Field(
        ..., description="The expression for this bond types group"
    )

    children: List[Union[ParametersUnitDef, BondType]] = Field(
        ..., description="Children of this bond type tag", alias="children"
    )

    def to_gmso_potentials(self, default_units):
        potentials = {"bond_types": {}}
        parameters_units = filter(
            lambda c: isinstance(c, ParametersUnitDef), self.children
        )
        units = {
            parameter_unit.parameter: u.Unit(parameter_unit.unit)
            for parameter_unit in parameters_units
        }

        for bond_type in filter(
            lambda c: isinstance(c, BondType), self.children
        ):
            bond_type_dict = bond_type.dict(
                by_alias=True,
                exclude={"children", "type1", "type2", "class1", "class2"},
                exclude_none=True,
            )

            if "expression" not in bond_type_dict:
                bond_type_dict["expression"] = self.expression

            if bond_type.type1 and bond_type.type2:
                bond_type_dict["member_types"] = (
                    bond_type.type1,
                    bond_type.type2,
                )

            elif bond_type.class1 and bond_type.class2:
                bond_type_dict["member_classes"] = (
                    bond_type.class1,
                    bond_type.class2,
                )

            bond_type_dict["parameters"] = bond_type.parameters(units)
            bond_type_dict["independent_variables"] = indep_vars(
                bond_type_dict["expression"],
                frozenset(bond_type_dict["parameters"]),
            )

            gmso_bond_type = GMSOBondType(**bond_type_dict)
            if gmso_bond_type.member_types:
                potentials["bond_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_bond_type.member_types)
                ] = gmso_bond_type
            else:
                potentials["bond_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_bond_type.member_classes)
                ] = gmso_bond_type

        return potentials

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        identifiers_registry = set()
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "BondType":
                bond_type = BondType.load_from_etree(el)
                identifier = tuple(
                    [bond_type.class1, bond_type.class2]
                    if bond_type.class1
                    else [bond_type.type1, bond_type.type2]
                )
                if identifier in identifiers_registry:
                    raise ValueError(
                        f"Duplicate entries found for BondType with identifiers {identifier}"
                    )
                else:
                    identifiers_registry.add(identifier)
                    identifiers_registry.add(tuple(reversed(identifier)))
                    children.append(bond_type)

        return cls(children=children, **attribs)


class AngleType(GMSOXMLTag):
    name: str = Field(
        None, description="The name of the angle type", alias="name"
    )

    class1: Optional[str] = Field(
        None, description="Class 1 for this angle type", alias="class1"
    )

    class2: Optional[str] = Field(
        None, description="Class 2 for this angle type", alias="class2"
    )

    class3: Optional[str] = Field(
        None, description="Class 3 for this angle type", alias="class3"
    )

    type1: Optional[str] = Field(
        None, description="Type 1 for this angle type", alias="type1"
    )

    type2: Optional[str] = Field(
        None, description="Type 2 for this angle type", alias="type2"
    )

    type3: Optional[str] = Field(
        None, description="Type 3 for this angle type", alias="type3"
    )

    children: List[Parameters] = Field(
        ..., description="The parameters and their values", alias="children"
    )

    @classmethod
    def load_from_etree(cls, root):
        children = []
        attribs = pad_with_wildcards(root.attrib, 3)
        for el in root.iterchildren():
            if el.tag == "Parameters":
                children.append(Parameters.load_from_etree(el))
        return cls(children=children, **attribs)


class AngleTypes(GMSOXMLChild):
    name: Optional[str] = Field(
        None, description="The name for this angle types group", alias="name"
    )

    expression: str = Field(
        ..., description="The expression for this angle types group"
    )

    children: List[Union[ParametersUnitDef, AngleType]] = Field(
        ..., description="Children of this angle types tag", alias="children"
    )

    def to_gmso_potentials(self, default_units):
        potentials = {"angle_types": {}}
        parameters_units = filter(
            lambda c: isinstance(c, ParametersUnitDef), self.children
        )
        units = {
            parameter_unit.parameter: u.Unit(parameter_unit.unit)
            for parameter_unit in parameters_units
        }

        for angle_type in filter(
            lambda c: isinstance(c, AngleType), self.children
        ):
            angle_type_dict = angle_type.dict(
                by_alias=True,
                exclude={
                    "children",
                    "type1",
                    "type2",
                    "type3",
                    "class1",
                    "class2",
                    "class3",
                },
                exclude_none=True,
            )

            if "expression" not in angle_type_dict:
                angle_type_dict["expression"] = self.expression

            if angle_type.type1 and angle_type.type2 and angle_type.type3:
                angle_type_dict["member_types"] = (
                    angle_type.type1,
                    angle_type.type2,
                    angle_type.type3,
                )

            elif angle_type.class1 and angle_type.class2 and angle_type.class3:
                angle_type_dict["member_classes"] = (
                    angle_type.class1,
                    angle_type.class2,
                    angle_type.class3,
                )

            angle_type_dict["parameters"] = angle_type.parameters(units)
            angle_type_dict["independent_variables"] = indep_vars(
                angle_type_dict["expression"],
                frozenset(angle_type_dict["parameters"]),
            )
            gmso_angle_type = GMSOAngleType(**angle_type_dict)
            if gmso_angle_type.member_types:
                potentials["angle_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_angle_type.member_types)
                ] = gmso_angle_type
            else:
                potentials["angle_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_angle_type.member_classes)
                ] = gmso_angle_type

        return potentials

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "AngleType":
                children.append(AngleType.load_from_etree(el))
        return cls(children=children, **attribs)


class TorsionType(GMSOXMLTag):
    name: str = Field(
        None, description="The name of the Dihedral/Improper type", alias="name"
    )

    class1: Optional[str] = Field(
        None,
        description="Class 1 for this Dihedral/Improper type",
        alias="class1",
    )

    class2: Optional[str] = Field(
        None,
        description="Class 2 for this Dihedral/Improper type",
        alias="class2",
    )

    class3: Optional[str] = Field(
        None,
        description="Class 3 for this Dihedral/Improper type",
        alias="class3",
    )

    class4: Optional[str] = Field(
        None,
        description="Class 4 for this Dihedral/Improper type",
        alias="class4",
    )

    type1: Optional[str] = Field(
        None,
        description="Type 1 for this Dihedral/Improper type",
        alias="type1",
    )

    type2: Optional[str] = Field(
        None,
        description="Type 2 for this Dihedral/Improper type",
        alias="type2",
    )

    type3: Optional[str] = Field(
        None,
        description="Type 3 for this Dihedral/Improper type",
        alias="type3",
    )

    type4: Optional[str] = Field(
        None,
        description="Type 4 for this Dihedral/Improper type",
        alias="type4",
    )

    children: List[Parameters] = Field(
        ..., description="The parameters and their values", alias="children"
    )

    @classmethod
    def load_from_etree(cls, root):
        children = []
        attribs = pad_with_wildcards(root.attrib, 4)
        for el in root.iterchildren():
            if el.tag == "Parameters":
                children.append(Parameters.load_from_etree(el))
        return cls(children=children, **attribs)


class DihedralType(TorsionType):
    pass


class ImproperType(TorsionType):
    pass


class TorsionTypes(GMSOXMLChild):
    name: Optional[str] = Field(
        None, description="The name for this angle types group", alias="name"
    )

    expression: str = Field(
        ..., description="The expression for this angle types group"
    )

    children: List[Union[ParametersUnitDef, TorsionType]] = Field(
        ...,
        description="Children of this dihedral/improper types tag",
        alias="children",
    )

    def to_gmso_potentials(self, default_units):
        potentials = {"dihedral_types": {}, "improper_types": {}}
        parameters_units = filter(
            lambda c: isinstance(c, ParametersUnitDef), self.children
        )
        units = {
            parameter_unit.parameter: u.Unit(parameter_unit.unit)
            for parameter_unit in parameters_units
        }

        for torsion_type in filter(
            lambda c: isinstance(c, (DihedralType, ImproperType)), self.children
        ):
            torsion_dict = torsion_type.dict(
                by_alias=True,
                exclude={
                    "children",
                    "type1",
                    "type2",
                    "type3",
                    "type4",
                    "class1",
                    "class2",
                    "class3",
                    "class4",
                },
                exclude_none=True,
            )

            if "expression" not in torsion_dict:
                torsion_dict["expression"] = self.expression

            if (
                torsion_type.type1
                and torsion_type.type2
                and torsion_type.type3
                and torsion_type.type4
            ):
                torsion_dict["member_types"] = (
                    torsion_type.type1,
                    torsion_type.type2,
                    torsion_type.type3,
                    torsion_type.type4,
                )

            elif (
                torsion_type.class1
                and torsion_type.class2
                and torsion_type.class3
                and torsion_type.class4
            ):
                torsion_dict["member_classes"] = (
                    torsion_type.class1,
                    torsion_type.class2,
                    torsion_type.class3,
                    torsion_type.class4,
                )

            torsion_dict["parameters"] = torsion_type.parameters(units)
            torsion_dict["independent_variables"] = indep_vars(
                torsion_dict["expression"],
                frozenset(torsion_dict["parameters"]),
            )
            if isinstance(torsion_type, DihedralType):
                gmso_torsion_type = GMSODihedralType(**torsion_dict)
                key = "dihedral_types"
            else:
                gmso_torsion_type = GMSOImproperType(**torsion_dict)
                key = "improper_types"

            if gmso_torsion_type.member_types:
                potentials[key][
                    FF_TOKENS_SEPARATOR.join(gmso_torsion_type.member_types)
                ] = gmso_torsion_type
            else:
                potentials[key][
                    FF_TOKENS_SEPARATOR.join(gmso_torsion_type.member_classes)
                ] = gmso_torsion_type

        return potentials

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "DihedralType":
                children.append(DihedralType.load_from_etree(el))
            elif el.tag == "ImproperType":
                child = ImproperType.load_from_etree(el)
                children.append(child)

        return cls(children=children, **attribs)


class ImproperTypes(TorsionTypes):
    pass


class DihedralTypes(TorsionTypes):
    pass


class PairPotentialType(GMSOXMLTag):
    name: str = Field(
        ..., description="Name of this PairPotential Type", alias="name"
    )

    type1: Optional[str] = Field(
        None, description="The type1 of this PairPotential Type", alias="type1"
    )

    type2: Optional[str] = Field(
        None, description="The type2 of this PairPotential Type", alias="type2"
    )

    class1: Optional[str] = Field(
        None,
        description="The class1 of this PairPotential Type",
        alias="class1",
    )

    class2: Optional[str] = Field(
        None,
        description="The class2 of this PairPotential Type",
        alias="class2",
    )

    children: List[Parameters] = Field(
        ..., description="The parameters and their values", alias="children"
    )

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            children.append(Parameters.load_from_etree(el))
        return cls(children=children, **attribs)


class PairPotentialTypes(GMSOXMLChild):
    name: Optional[str] = Field(
        None,
        description="The name of this pair potential types group",
        alias="name",
    )

    expression: str = Field(
        ...,
        description="The expression for this pair potential types group",
        alias="expression",
    )

    children: List[Union[ParametersUnitDef, PairPotentialType]] = Field(
        ..., description="The children", alias="children"
    )

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "PairPotentialType":
                children.append(PairPotentialType.load_from_etree(el))
        return cls(children=children, **attribs)


class Units(GMSOXMLTag):
    energy: Optional[str] = Field(None, alias="energy")

    distance: Optional[str] = Field(None, alias="distance")

    mass: Optional[str] = Field(None, alias="mass")

    charge: Optional[str] = Field(None, alias="charge")

    temperature: Optional[str] = Field(None, alias="temperature")

    angle: Optional[str] = Field(None, alias="angle")

    time: Optional[str] = Field(None, alias="time")

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        return cls(**attribs)


class FFMetaData(GMSOXMLChild):
    children: List[Units] = Field([], alias="children")

    electrostatics14Scale: float = Field(0.5, alias="electrostatics14Scale")

    nonBonded14Scale: float = Field(0.5, alias="nonBonded14Scale")

    combining_rule: str = Field("geometric", alias="combiningRule")

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for unit in root.iterchildren():
            children.append(Units.load_from_etree(unit))
        return cls(children=children, **attribs)

    def gmso_scaling_factors(self):
        return self.dict(
            include={"electrostatics14Scale", "nonBonded14Scale"},
            exclude_none=True,
        )

    def get_default_units(self):
        units_dict = {}
        units = self.children[0].dict(by_alias=True, exclude_none=True)
        for name, unit in units.items():
            try:
                units_dict[name] = u.Unit(unit)
            except u.exceptions.UnitParseError:
                units_dict[name] = getattr(u.physical_constants, unit)

        return units_dict


class ForceField(GMSOXMLTag):
    name: str = Field(
        "ForceField", description="Name of the ForceField", alias="name"
    )

    version: str = Field(
        "1.0.0", description="The version of the ForceField", alias="version"
    )

    children: List[GMSOXMLChild] = Field(
        ..., description="The children tags", alias="children"
    )

    def to_gmso_ff(self):
        ff = GMSOForceField()
        ff.name = self.name
        ff.version = self.version
        metadata = list(
            filter(lambda child: isinstance(child, FFMetaData), self.children)
        ).pop()
        default_units = metadata.get_default_units()
        ff.scaling_factors = metadata.gmso_scaling_factors()
        ff.combining_rule = metadata.combining_rule
        remaining_children = filter(
            lambda c: not isinstance(c, (FFMetaData, Units)),
            self.children,
        )
        ff_potentials = {}

        for child in remaining_children:
            if hasattr(child, "to_gmso_potentials"):
                potentials = child.to_gmso_potentials(default_units)
                for attr in potentials:
                    if attr in ff_potentials:
                        ff_potentials[attr].update(potentials[attr])
                    else:
                        ff_potentials[attr] = potentials[attr]

        for attr in ff_potentials:
            setattr(ff, attr, ff_potentials[attr])

        return ff

    @classmethod
    def load_from_etree(cls, root) -> "ForceField":
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "FFMetaData":
                children.append(FFMetaData.load_from_etree(el))
            if el.tag == "AtomTypes":
                children.append(AtomTypes.load_from_etree(el))
            elif el.tag == "BondTypes":
                children.append(BondTypes.load_from_etree(el))
            elif el.tag == "AngleTypes":
                children.append(AngleTypes.load_from_etree(el))
            elif el.tag == "DihedralTypes":
                children.append(DihedralTypes.load_from_etree(el))
            elif el.tag == "ImproperTypes":
                children.append(ImproperTypes.load_from_etree(el))
            elif el.tag == "PairPotentialTypes":
                children.append(PairPotentialTypes.load_from_etree(el))
        return cls(children=children, **attribs)
