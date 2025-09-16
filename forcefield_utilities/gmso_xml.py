import re
from functools import lru_cache
from typing import List, Optional, Set, Tuple, Union

import numpy as np
import sympy
import unyt as u
from gmso import ForceField as GMSOForceField
from gmso.core.angle_type import AngleType as GMSOAngleType
from gmso.core.atom_type import AtomType as GMSOAtomType
from gmso.core.bond_type import BondType as GMSOBondType
from gmso.core.dihedral_type import DihedralType as GMSODihedralType
from gmso.core.improper_type import ImproperType as GMSOImproperType
from gmso.core.pairpotential_type import (
    PairPotentialType as GMSOPairPotentialType,
)
from gmso.core.virtual_type import (
    VirtualPositionType as GMSOVirtualPositionType,
)
from gmso.core.virtual_type import (
    VirtualPotentialType as GMSOVirtualPotentialType,
)
from gmso.core.virtual_type import VirtualType as GMSOVirtualType
from gmso.utils._constants import FF_TOKENS_SEPARATOR
from gmso.utils.ff_utils import _get_member_classes, _get_member_types
from pydantic import BaseModel, ConfigDict, Field

# TODO: add custom unyt registry
from unyt import Unit, UnitRegistry

from forcefield_utilities.utils import (
    get_virtual_ntype_or_nclass_attribs,
    pad_with_wildcards,
)

reg = UnitRegistry()
charge_dim = u.dimensions.current_mks * u.dimensions.time
elementary_charge_conversion = (
    1 * getattr(u.physical_constants, "elementary_charge").value
)
reg.add(
    "elementary_charge",
    base_value=elementary_charge_conversion,
    dimensions=charge_dim,
    tex_repr=r"\rm{e}",
)

kb_dim = u.dimensions.energy / u.dimensions.temperature
kb_conversion = (
    1 * getattr(u.physical_constants, "boltzmann_constant_mks").value
)
reg.add("kb", base_value=kb_conversion, dimensions=kb_dim, tex_repr=r"\rm{kb}")


def get_identifiers_registry():
    return {
        "AtomTypes": set(),
        "BondTypes": set(),
        "AngleTypes": set(),
        "DihedralTypes": set(),
        "ImproperTypes": set(),
        "PairPotentialTypes": set(),
    }


def reverse_identifier(identifier: str):
    bond_tokens = ["~", "-", "=", "#"]
    outStr = ""
    currentNode = ""
    for letter in identifier[::-1]:
        if letter in bond_tokens:
            outStr += currentNode[::-1]
            currentNode = ""
            outStr += letter  # should be a bond
        else:
            currentNode += letter
    if currentNode:
        outStr += currentNode[::-1]
    return outStr


def register_identifiers(registry, identifier, for_type="AtomTypes"):
    if identifier in registry:
        raise ValueError(
            f"Duplicate identifier found for {for_type}: {identifier}"
        )

    if for_type == "AtomTypes":
        registry.add(identifier)
    elif (
        for_type == "BondTypes"
        or for_type == "AngleTypes"
        or for_type == "DihedralTypes"
        or for_type == "PairPotentialTypes"
    ):
        if isinstance(identifier, str):
            registry.add(identifier)
            registry.add(reverse_identifier(identifier))
        else:
            identifierStr = "".join(identifier)
            registry.add(identifierStr)
            registry.add(reverse_identifier(identifierStr))
    elif for_type == "ImproperTypes":
        if isinstance(identifier, str):
            (central, second, third, fourth) = re.split(
                "(?=[\~\-\=\#])", identifier
            )
        else:  # identifier is a tuple
            (central, second, third, fourth) = identifier
        mirrors = [
            (central, second, third, fourth),
            (central, second, fourth, third),
            (central, third, second, fourth),
            (central, third, fourth, second),
            (central, fourth, second, third),
            (central, fourth, third, second),
        ]
        for mirror in mirrors:
            mirrorStr = "".join(mirror)
            registry.add(mirrorStr)


@lru_cache(maxsize=128)
def indep_vars(expr: str, dependent: frozenset) -> Set:
    """Given an expression and dependent variables, return independent variables for it."""
    dependent_symbols = frozenset(map(lambda x: sympy.symbols(x), dependent))
    return sympy.sympify(expr).free_symbols - dependent_symbols


class GMSOXMLTag(BaseModel):
    """The Base GMSO XML Class. Used for convience."""

    model_config = ConfigDict(
        populate_by_name=True,
        frozen=True,
        arbitrary_types_allowed=True,
    )

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

    value: Union[np.ndarray, float] = Field(
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
            parameter_unit.parameter: u.Unit(parameter_unit.unit, registry=reg)
            for parameter_unit in parameters_units
        }

        for atom_type in filter(
            lambda c: isinstance(c, AtomType), self.children
        ):
            atom_type_dict = atom_type.model_dump(
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
    def load_from_etree(cls, root, existing):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "AtomType":
                atom_type = AtomType.load_from_etree(el)
                identifier = atom_type.name
                register_identifiers(existing, identifier, "AtomTypes")
                children.append(atom_type)
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

    classes: Optional[str] = Field(
        None,
        description="Identifying string for both classes in this bond type, with bond order specified as ~,-,=,#",
        alias="classes",
    )

    type1: Optional[str] = Field(
        None, description="Type 1 for this bond type", alias="type1"
    )

    type2: Optional[str] = Field(
        None, description="Type 2 for this bond type", alias="type2"
    )

    types: Optional[str] = Field(
        None,
        description="Identifying string for both types in this bond type, with bond order specified as ~,-,=,#",
        alias="types",
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
            bond_type_dict = bond_type.model_dump(
                by_alias=True,
                exclude={
                    "children",
                    "type1",
                    "type2",
                    "class1",
                    "class2",
                    "types",
                    "classes",
                },
                exclude_none=True,
            )

            if "expression" not in bond_type_dict:
                bond_type_dict["expression"] = self.expression

            identifier = None
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
            elif bond_type.types:
                type1, type2 = re.split(r"[\~\-\=\#]+", bond_type.types)
                bond_type_dict["member_types"] = (
                    type1,
                    type2,
                )
                identifier = bond_type.types

            elif bond_type.classes:
                class1, class2 = re.split(r"[\~\-\=\#]+", bond_type.classes)
                bond_type_dict["member_classes"] = (
                    class1,
                    class2,
                )
                identifier = bond_type.classes

            bond_type_dict["parameters"] = bond_type.parameters(units)
            bond_type_dict["independent_variables"] = indep_vars(
                bond_type_dict["expression"],
                frozenset(bond_type_dict["parameters"]),
            )

            gmso_bond_type = GMSOBondType(**bond_type_dict)
            if identifier:
                potentials["bond_types"][identifier] = gmso_bond_type

            elif gmso_bond_type.member_types:  # create wildcard identifier
                potentials["bond_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_bond_type.member_types)
                ] = gmso_bond_type
            else:
                potentials["bond_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_bond_type.member_classes)
                ] = gmso_bond_type

        return potentials

    @classmethod
    def load_from_etree(cls, root, existing):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "BondType":
                bond_type = BondType.load_from_etree(el)
                if bond_type.types:
                    identifier = bond_type.types
                elif bond_type.classes:
                    identifier = bond_type.classes
                elif bond_type.type1:
                    identifier = f"{bond_type.type1}~"
                elif bond_type.class1:
                    identifier = f"{bond_type.class1}~"
                if bond_type.type2:
                    identifier += bond_type.type2
                elif bond_type.class2:
                    identifier += bond_type.class2
                register_identifiers(existing, identifier, "BondTypes")
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

    classes: Optional[str] = Field(
        None,
        description="Identifying string for all classes in this angle type, with bond order specified as ~,-,=,#",
        alias="classes",
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

    types: Optional[str] = Field(
        None,
        description="Identifying string for all types in this angle type, with bond order specified as ~,-,=,#",
        alias="types",
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
            angle_type_dict = angle_type.model_dump(
                by_alias=True,
                exclude={
                    "children",
                    "type1",
                    "type2",
                    "type3",
                    "class1",
                    "class2",
                    "class3",
                    "types",
                    "classes",
                },
                exclude_none=True,
            )

            if "expression" not in angle_type_dict:
                angle_type_dict["expression"] = self.expression

            identifier = None
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
            elif angle_type.types:
                # TODO: Set bond orders
                type1, type2, type3 = re.split(r"[\~\-\=\#]+", angle_type.types)
                angle_type_dict["member_types"] = (
                    type1,
                    type2,
                    type3,
                )
                identifier = angle_type.types

            elif angle_type.classes:
                # TODO: Set bond orders
                class1, class2, class3 = re.split(
                    r"[\~\-\=\#]+", angle_type.classes
                )
                angle_type_dict["member_classes"] = (
                    class1,
                    class2,
                    class3,
                )
                identifier = angle_type.classes

            angle_type_dict["parameters"] = angle_type.parameters(units)
            angle_type_dict["independent_variables"] = indep_vars(
                angle_type_dict["expression"],
                frozenset(angle_type_dict["parameters"]),
            )
            gmso_angle_type = GMSOAngleType(**angle_type_dict)
            if identifier:
                potentials["angle_types"][identifier] = gmso_angle_type
            elif gmso_angle_type.member_types:
                potentials["angle_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_angle_type.member_types)
                ] = gmso_angle_type
            else:
                potentials["angle_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_angle_type.member_classes)
                ] = gmso_angle_type

        return potentials

    @classmethod
    def load_from_etree(cls, root, existing):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "AngleType":
                angle_type = AngleType.load_from_etree(el)
                if angle_type.types:
                    identifier = angle_type.types
                elif angle_type.classes:
                    identifier = angle_type.classes
                elif angle_type.type1:
                    identifier = f"{angle_type.type1}~"
                elif angle_type.class1:
                    identifier = f"{angle_type.class1}~"
                if angle_type.type2:
                    identifier += angle_type.type2 + "~"
                elif angle_type.class2:
                    identifier += angle_type.class2 + "~"
                if angle_type.type3:
                    identifier += angle_type.type3
                elif angle_type.class3:
                    identifier += angle_type.class3
                register_identifiers(existing, identifier, "AngleTypes")
                children.append(angle_type)
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

    classes: Optional[str] = Field(
        None,
        description="Identifying string for all classes in this torsion type, with bond order specified as ~,-,=,#",
        alias="classes",
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

    types: Optional[str] = Field(
        None,
        description="Identifying string for all types in this torsion type, with bond order specified as ~,-,=,#",
        alias="types",
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
            torsion_dict = torsion_type.model_dump(
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
                    "types",
                    "classes",
                },
                exclude_none=True,
            )

            if "expression" not in torsion_dict:
                torsion_dict["expression"] = self.expression

            identifier = None
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
            elif torsion_type.types:
                # TODO: Set bond orders
                type1, type2, type3, type4 = re.split(
                    r"[\~\-\=\#]+", torsion_type.types
                )
                torsion_dict["member_types"] = (
                    type1,
                    type2,
                    type3,
                    type4,
                )
                identifier = torsion_type.types

            elif torsion_type.classes:
                # TODO: Set bond orders
                class1, class2, class3, class4 = re.split(
                    r"[\~\-\=\#]+", torsion_type.classes
                )
                torsion_dict["member_classes"] = (
                    class1,
                    class2,
                    class3,
                    class4,
                )
                identifier = torsion_type.classes

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

            if identifier:
                potentials[key][identifier] = gmso_torsion_type
            elif gmso_torsion_type.member_types:
                potentials[key][
                    FF_TOKENS_SEPARATOR.join(gmso_torsion_type.member_types)
                ] = gmso_torsion_type
            else:
                potentials[key][
                    FF_TOKENS_SEPARATOR.join(gmso_torsion_type.member_classes)
                ] = gmso_torsion_type

        return potentials

    @classmethod
    def load_from_etree(cls, root, existing_dihedrals, existing_impropers):
        attribs = root.attrib
        children = []
        child_loaders = {
            "DihedralType": DihedralType,
            "ImproperType": ImproperType,
        }
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "DihedralType" or el.tag == "ImproperType":
                tor_type = child_loaders[el.tag].load_from_etree(el)
                if tor_type.types:
                    identifier = tor_type.types
                elif tor_type.classes:
                    identifier = tor_type.classes
                elif tor_type.type1:
                    identifier = f"{tor_type.type1}~"
                elif tor_type.class1:
                    identifier = f"{tor_type.class1}~"
                if tor_type.type2:
                    identifier += tor_type.type2 + "~"
                elif tor_type.class2:
                    identifier += tor_type.class2 + "~"
                if tor_type.type3:
                    identifier += tor_type.type3 + "~"
                elif tor_type.class3:
                    identifier += tor_type.class3 + "~"
                if tor_type.type4:
                    identifier += tor_type.type4
                elif tor_type.class4:
                    identifier += tor_type.class4
                register_identifiers(
                    (
                        existing_impropers
                        if el.tag == "ImproperType"
                        else existing_dihedrals
                    ),
                    identifier,
                    el.tag + "s",
                )
                children.append(tor_type)

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
        attribs = pad_with_wildcards(root.attrib, 2)
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
    def load_from_etree(cls, root, existing):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "PairPotentialType":
                pptype = PairPotentialType.load_from_etree(el)
                identifier = tuple(
                    [pptype.class1, pptype.class2]
                    if pptype.class1
                    else [pptype.type1, pptype.type2]
                )
                register_identifiers(existing, identifier, "PairPotentialTypes")
                children.append(pptype)
        return cls(children=children, **attribs)

    def to_gmso_potentials(self, default_units):
        potentials = {"pairpotential_types": {}}
        parameters_units = filter(
            lambda c: isinstance(c, ParametersUnitDef), self.children
        )
        units = {
            parameter_unit.parameter: u.Unit(parameter_unit.unit)
            for parameter_unit in parameters_units
        }

        for pairpotential_type in filter(
            lambda c: isinstance(c, PairPotentialType), self.children
        ):
            pairpotential_type_dict = pairpotential_type.model_dump(
                by_alias=True,
                exclude={"children", "type1", "type2", "class1", "class2"},
                exclude_none=True,
            )

            if "expression" not in pairpotential_type_dict:
                pairpotential_type_dict["expression"] = self.expression

            if pairpotential_type.type1 and pairpotential_type.type2:
                pairpotential_type_dict["member_types"] = (
                    pairpotential_type.type1,
                    pairpotential_type.type2,
                )

            elif pairpotential_type.class1 and pairpotential_type.class2:
                pairpotential_type_dict["member_classes"] = (
                    pairpotential_type.class1,
                    pairpotential_type.class2,
                )

            pairpotential_type_dict["parameters"] = (
                pairpotential_type.parameters(units)
            )
            pairpotential_type_dict["independent_variables"] = indep_vars(
                pairpotential_type_dict["expression"],
                frozenset(pairpotential_type_dict["parameters"]),
            )

            gmso_pairpotential_type = GMSOPairPotentialType(
                **pairpotential_type_dict
            )
            if gmso_pairpotential_type.member_types:
                potentials["pairpotential_types"][
                    FF_TOKENS_SEPARATOR.join(
                        gmso_pairpotential_type.member_types
                    )
                ] = gmso_pairpotential_type
            else:
                potentials["pairpotential_types"][
                    FF_TOKENS_SEPARATOR.join(
                        gmso_pairpotential_type.member_classes
                    )
                ] = gmso_pairpotential_type

        return potentials


class VirtualPotentialType(GMSOXMLTag):
    children: List[Parameters] = Field(
        ...,
        description="The unyt parameters and their values",
        alias="children",
    )

    @classmethod
    def load_from_etree(cls, root):
        attribs = pad_with_wildcards(root.attrib, 2)
        children = []
        for el in root.iterchildren():
            children.append(Parameters.load_from_etree(el))
        return cls(children=children, **attribs)


class VirtualPositionType(GMSOXMLTag):
    children: List[Parameters] = Field(
        ...,
        description="The unyt parameters and their values",
        alias="children",
    )

    @classmethod
    def load_from_etree(cls, root):
        attribs = pad_with_wildcards(root.attrib, 2)
        children = []
        for el in root.iterchildren():
            children.append(Parameters.load_from_etree(el))
        return cls(children=children, **attribs)


class VirtualSiteType(GMSOXMLTag):
    name: str = Field(
        None, description="The name of the virtual site type", alias="name"
    )

    member_types: Optional[Tuple[str, ...]] = Field(
        None,
        description="List-like of gmso.AtomType.name "
        "defining the members of this angle type",
        alias="member_types",
        min_length=0,
        max_length=12,
    )

    member_classes: Optional[Tuple[str, ...]] = Field(
        None,
        description="List-like of gmso.AtomType.atomclass "
        "defining the members of this angle type",
        alias="member_classes",
        min_length=0,
        max_length=12,
    )

    charge: Optional[float] = Field(
        None, description="The charge of the virtual type", alias="charge"
    )

    virtual_potential: Optional[VirtualPotentialType] = Field(
        default=None,
        description="virtual type defining the interaction energy for a virtual site",
        alias="virtual_potential",
    )

    virtual_position: Optional[VirtualPositionType] = Field(
        default=None,
        description="virtual type defining the position for a virtual site",
        alias="virtual_position",
    )

    doi: Optional[str] = Field(
        None, description="The doi of this virtual type", alias="doi"
    )

    @classmethod
    def load_from_etree(cls, root):
        attribs = get_virtual_ntype_or_nclass_attribs(root.attrib)
        for el in root.iterchildren():
            if el.tag == "Potential":
                attribs["virtual_potential"] = (
                    VirtualPotentialType.load_from_etree(el)
                )
            elif el.tag == "Position":
                attribs["virtual_position"] = (
                    VirtualPositionType.load_from_etree(el)
                )
        return cls(**attribs)


class VirtualSiteTypes(GMSOXMLChild):
    name: Optional[str] = Field(
        None, description="The name for this virtual types group", alias="name"
    )

    potential_expression: str = Field(
        ..., description="The general expression for the potential energy"
    )

    position_expression: str = Field(
        ...,
        description="The general expression for calculating the position of the virtual site from its parent sites",
    )

    children: List[Union[VirtualSiteType, ParametersUnitDef]] = Field(
        ..., description="Children of this virtual types tag", alias="children"
    )

    def to_gmso_potentials(self, default_units):
        potentials = {"virtual_types": {}}
        parameters_units = filter(
            lambda c: isinstance(c, ParametersUnitDef), self.children
        )
        units = {
            parameter_unit.parameter: u.Unit(parameter_unit.unit)
            for parameter_unit in parameters_units
        }

        for virtual_type in filter(
            lambda c: isinstance(c, VirtualSiteType), self.children
        ):
            virtual_type_dict = virtual_type.model_dump(
                by_alias=True,
                exclude={
                    "children",
                    "member_types",
                    "member_classes",
                    "virtual_potential",
                    "virtual_position",
                },
                exclude_none=True,
            )
            if virtual_type.member_types:
                virtual_type_dict["member_types"] = virtual_type.member_types

            elif virtual_type.member_classes:
                virtual_type_dict["member_classes"] = (
                    virtual_type.member_classes
                )

            potentialDict = {}
            if self.potential_expression:
                potentialDict["expression"] = self.potential_expression
            if virtual_type.virtual_potential:
                potentialDict["parameters"] = (
                    virtual_type.virtual_potential.parameters(units)
                )
            potentialDict["independent_variables"] = indep_vars(
                potentialDict["expression"],
                frozenset(potentialDict["parameters"]),
            )
            virtual_type_dict["virtual_potential"] = GMSOVirtualPotentialType(
                **potentialDict
            )

            positionDict = {}
            if self.position_expression:
                positionDict["expression"] = self.position_expression
            if virtual_type.virtual_position:
                positionDict["parameters"] = (
                    virtual_type.virtual_position.parameters(units)
                )
            positionDict["independent_variables"] = indep_vars(
                positionDict["expression"],
                frozenset(positionDict["parameters"]),
            )
            virtual_type_dict["virtual_position"] = GMSOVirtualPositionType(
                **positionDict
            )
            gmso_virtual_type = GMSOVirtualType(**virtual_type_dict)
            if gmso_virtual_type.member_types:
                potentials["virtual_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_virtual_type.member_types)
                ] = gmso_virtual_type
            else:
                potentials["virtual_types"][
                    FF_TOKENS_SEPARATOR.join(gmso_virtual_type.member_classes)
                ] = gmso_virtual_type

        return potentials

    @classmethod
    def load_from_etree(cls, root, existing):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "Potential":
                for child in el.iterfind("ParametersUnitDef"):
                    children.append(ParametersUnitDef.load_from_etree(child))
                attribs.update({"potential_expression": el.get("expression")})
            elif el.tag == "Position":
                for child in el.iterfind("ParametersUnitDef"):
                    children.append(ParametersUnitDef.load_from_etree(child))
                attribs.update({"position_expression": el.get("expression")})
            elif el.tag == "VirtualSiteType":
                virtual_type = VirtualSiteType.load_from_etree(el)
                member_types = _get_member_types(el)
                member_classes = _get_member_classes(el)
                identifier = tuple(
                    member_types if member_types else member_classes
                )
                register_identifiers(existing, identifier, "VirtualSiteTypes")
                children.append(virtual_type)
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
        return self.model_dump(
            include={"electrostatics14Scale", "nonBonded14Scale"},
            exclude_none=True,
        )

    def get_default_units(self):
        units_dict = {}
        units = self.children[0].model_dump(by_alias=True, exclude_none=True)
        for name, unit in units.items():
            unit_object = source_units(unit)
            if unit_object:
                units_dict[name] = unit_object
            else:
                raise u.exceptions.UnitParseError(
                    f"Unit {unit} with name {name} in the forcefield {self.name} cannot be handled. \
                    Please consider adding it to the unit registry."
                )

        return units_dict


def source_units(unit):
    try:
        attach_unit = u.Unit(unit, registry=reg)
    except u.exceptions.UnitParseError:
        try:
            attach_unit = Unit("electron_charge", registry=reg)
        except u.exceptions.UnitParseError:
            attach_unit = getattr(u.physical_constants, unit)
    return attach_unit


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
        ff.units = default_units
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
        identifiers_registry = get_identifiers_registry()
        for el in root.iterchildren():
            if el.tag == "FFMetaData":
                children.append(FFMetaData.load_from_etree(el))
            if el.tag == "AtomTypes":
                children.append(
                    AtomTypes.load_from_etree(
                        el, identifiers_registry["AtomTypes"]
                    )
                )
            elif el.tag == "BondTypes":
                children.append(
                    BondTypes.load_from_etree(
                        el, identifiers_registry["BondTypes"]
                    )
                )
            elif el.tag == "AngleTypes":
                children.append(
                    AngleTypes.load_from_etree(
                        el, identifiers_registry["AngleTypes"]
                    )
                )
            elif el.tag == "DihedralTypes":
                children.append(
                    DihedralTypes.load_from_etree(
                        el,
                        identifiers_registry["DihedralTypes"],
                        identifiers_registry["ImproperTypes"],
                    )
                )
            elif el.tag == "ImproperTypes":
                children.append(
                    ImproperTypes.load_from_etree(
                        el,
                        identifiers_registry["DihedralTypes"],
                        identifiers_registry["ImproperTypes"],
                    )
                )
            elif el.tag == "PairPotentialTypes":
                children.append(
                    PairPotentialTypes.load_from_etree(
                        el, identifiers_registry["PairPotentialTypes"]
                    )
                )
            elif el.tag == "VirtualSiteTypes":
                children.append(
                    VirtualSiteTypes.load_from_etree(
                        el, identifiers_registry["PairPotentialTypes"]
                    )
                )

        return cls(children=children, **attribs)
