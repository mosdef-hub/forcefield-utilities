from typing import ClassVar, List, Optional, Set

from lxml.etree import ElementTree
from pydantic import BaseModel, Field

__all__ = ["ForceField"]


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


class Bond(ForceFieldChild):
    class1: Optional[str] = Field(
        default=None, description="Class 1 of the bond", alias="class1"
    )

    class2: Optional[str] = Field(
        default=None, description="Class 2 of the bond", alias="class2"
    )

    type1: Optional[str] = Field(
        default=None, description="Type 1 of the bond", alias="class2"
    )

    type2: Optional[str] = Field(
        default=None, description="Type 2 of the bond", alias="type2"
    )

    length: Optional[float] = Field(
        default=None, description="The length of the bond", alias="length"
    )


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


class HarmonicBondForce(ForceFieldChild):
    pass


child_mappers = {"AtomTypes": AtomTypes, "HarmonicBondForce": HarmonicBondForce}


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
            if el.tag in child_mappers:
                children.append(child_mappers[el.tag].load_from_etree(el))
        return cls(children=children, **attribs)
