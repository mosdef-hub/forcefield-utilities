from functools import wraps
from typing import ClassVar, List, Optional, Set

from pydantic import BaseModel, Field

__all__ = ["ForceField"]

loaders = {}


class registers_loader:
    def __init__(self, name):
        self.name = name

    def __call__(self, child_class):
        loaders[self.name] = child_class


class FoyerXMLTag(BaseModel):
    """The Base Foyer XML Class. Used for convience"""

    def to_etree(self):
        pass

    class Config:
        allow_population_by_field_name = True
        allow_mutation = False
        arbitrary_types_allowed = True


class ForceFieldChild(FoyerXMLTag):
    """Any XML Child of <ForceField>"""

    children: Optional[Set["ForceFieldChild"]] = None


class FoyerXMLAtomType(ForceFieldChild):
    pass


class Type(FoyerXMLAtomType):
    name: str = Field(..., description="The AtomType Name", alias="name")

    atom_class: str = Field(
        ..., description="The AtomType class", alias="class"
    )

    element: Optional[str] = Field(
        default=None,
        description="The element of of the AtomType",
        alias="element",
    )

    mass: float = Field(
        ..., description="The mass of the AtomType", alias="mass"
    )

    smarts_def: Optional[str] = Field(
        default=None, description="The smarts definition", alias="def"
    )

    desc: Optional[str] = Field(
        default=None, description="Description of the atom type", alias="desc"
    )

    overrides: Optional[str] = Field(
        default=None, description="The overrides", alias="overrides"
    )

    doi: Optional[str] = Field(default=None, alias="doi", description="The doi")


@registers_loader(name="AtomTypes")
class AtomTypes(ForceFieldChild):
    children: List[Type] = Field(
        ..., description="The AtomType definitions", alias="types"
    )

    @classmethod
    def load_from_etree(cls, atom_types):
        children = []
        for atom_type in atom_types.iterchildren():
            if atom_type.tag == Type.__name__:
                children.append(Type.parse_obj(atom_type.attrib))
        return cls(children=children)


class Bond(ForceFieldChild):
    class1: Optional[str] = Field(
        default=None, description="Class 1 of the bond", alias="class1"
    )

    class2: Optional[str] = Field(
        default=None, description="Class 2 of the bond", alias="class2"
    )

    type1: Optional[str] = Field(
        default=None, description="Type 1 of the bond", alias="type1"
    )

    type2: Optional[str] = Field(
        default=None, description="Type 2 of the bond", alias="type2"
    )

    length: float = Field(
        ..., description="The length of the bond", alias="length"
    )

    k: float = Field(..., description="The k-value of the bond", alias="k")


@registers_loader(name="HarmonicBondForce")
class HarmonicBondForce(ForceFieldChild):
    children: List[Bond]

    @classmethod
    def load_from_etree(cls, bonds):
        children = []
        for bond_type in bonds.iterchildren():
            if bond_type.tag == Bond.__name__:
                children.append(Bond.parse_obj(bond_type.attrib))
        return cls(children=children)


class Angle(ForceFieldChild):
    class1: Optional[str] = Field(
        default=None, description="Class 1 of the angle", alias="class1"
    )

    class2: Optional[str] = Field(
        default=None, description="Class 2 of the angle", alias="class2"
    )

    class3: Optional[str] = Field(
        default=None, description="Class 3 of the angle", alias="class3"
    )

    type1: Optional[str] = Field(
        default=None, description="Type 1 of the angle", alias="type1"
    )

    type2: Optional[str] = Field(
        default=None, description="Type 2 of the angle", alias="type2"
    )

    type3: Optional[str] = Field(
        default=None, description="Type 3 of the angle", alias="type3"
    )

    angle: float = Field(..., description="The angle", alias="angle")

    k: float = Field(..., description="The k-value of angle", alias="k")


@registers_loader(name="HarmonicAngleForce")
class HarmonicAngleForce(ForceFieldChild):
    children: List[Angle] = []

    @classmethod
    def load_from_etree(cls, angles):
        children = []
        for angle_type in angles.iterchildren():
            if angle_type.tag == Angle.__name__:
                children.append(Angle.parse_obj(angle_type.attrib))
        return cls(children=children)


class ForceField(FoyerXMLTag):
    name: str = Field(
        default="Forcefield", alias="name", description="Name of the Forcefield"
    )

    version: str = Field(
        default="1.0.0",
        alias="version",
        description="Version of the Forcefield",
    )

    combining_rule: str = Field(
        default="geometric",
        alias="combining_rule",
        description="The Combining rule",
    )

    children: List[ForceFieldChild] = Field(
        ..., alias="children", description="The children tags"
    )

    def to_foyer_etree(self):
        pass

    @classmethod
    def load_from_etree(cls, root) -> "ForceField":
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag in loaders:
                children.append(loaders[el.tag].load_from_etree(el))
        return cls(children=children, **attribs)
