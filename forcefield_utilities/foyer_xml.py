import re
import warnings
from typing import Any, List, Optional, Set, Union

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


class Dihedral(ForceFieldChild):
    class1: Optional[str] = Field(
        default=None, description="Class 1 of the dihedral", alias="class1"
    )

    class2: Optional[str] = Field(
        default=None, description="Class 2 of the dihedral", alias="class2"
    )

    class3: Optional[str] = Field(
        default=None, description="Class 3 of the dihedral", alias="class3"
    )

    class4: Optional[str] = Field(
        default=None, description="Class 4 of the dihedral", alias="class4"
    )

    type1: Optional[str] = Field(
        default=None, description="Type 1 of the dihedral", alias="type1"
    )

    type2: Optional[str] = Field(
        default=None, description="Type 2 of the dihedral", alias="type2"
    )

    type3: Optional[str] = Field(
        default=None, description="Type 3 of the dihedral", alias="type3"
    )

    type4: Optional[str] = Field(
        default=None, description="Type 4 of the dihedral", alias="type4"
    )


class RBDihedral(Dihedral):
    c0: float = Field(..., description="C0 Parameter", alias="c0")

    c1: float = Field(..., description="C1 Parameter", alias="c1")

    c2: float = Field(..., description="C2 Parameter", alias="c2")

    c3: float = Field(..., description="C3 Parameter", alias="c3")

    c4: float = Field(..., description="C4 Parameter", alias="c4")

    c5: float = Field(..., description="C5 Parameter", alias="c5")


class RBProper(RBDihedral):
    pass


class RBImproper(RBDihedral):
    pass


@registers_loader(name="RBTorsionForce")
class RBTorsionForce(ForceFieldChild):
    children: List[Union[RBProper, RBImproper]]

    @classmethod
    def load_from_etree(cls, torsions):
        children = []
        for dihedral_type in torsions.iterchildren():
            if dihedral_type.tag == "Proper":
                Creator = RBProper
            elif dihedral_type.tag == "Improper":
                Creator = RBImproper
            else:
                warnings.warn(
                    f"Tag {dihedral_type.tag} not understood skipping"
                )
                continue
            children.append(Creator.parse_obj(dihedral_type.attrib))

        return cls(children=children)


class PeriodicDihedral(Dihedral):
    periodicity: List[float] = Field(
        ..., description="Periodicity 1, 2, ...", alias="periodicity"
    )

    k: List[float] = Field(..., description="Phase 1, 2, ....", alias="k")

    phase: List[float] = Field(
        ..., description="Phase 1, 2, 3 ...", alias="phase"
    )

    @staticmethod
    def periodic_attribs_to_list(attrib):
        max_periodicity = 1
        for key in attrib:
            if key.startswith("periodicity"):
                periodicity_count = filter(
                    lambda x: bool(x), re.split(r"[a-zA-z]+", key)
                )
                periodicity_count = int(list(periodicity_count).pop())
                if periodicity_count and periodicity_count > max_periodicity:
                    max_periodicity = periodicity_count
        attrib_dict = {"periodicity": [], "phase": [], "k": []}

        for j in range(1, max_periodicity + 1):
            attrib_dict["periodicity"].append(attrib[f"periodicity{j}"])
            attrib_dict["phase"].append(attrib[f"phase{j}"])
            attrib_dict["k"].append(attrib[f"k{j}"])
        return attrib_dict

    @classmethod
    def parse_obj(cls, obj: Any):
        dict_attribs = cls.periodic_attribs_to_list(obj)
        return super().parse_obj(dict_attribs)


class PeriodicProper(PeriodicDihedral):
    pass


class PeriodicImproper(PeriodicDihedral):
    pass


@registers_loader(name="PeriodicTorsionForce")
class PeriodicTorsionForce(ForceFieldChild):
    children: List[PeriodicDihedral]

    @classmethod
    def load_from_etree(cls, torsions):
        children = []
        for dihedral_type in torsions.iterchildren():
            if dihedral_type.tag == "Proper":
                Creator = PeriodicProper
            elif dihedral_type.tag == "Improper":
                Creator = PeriodicImproper
            else:
                warnings.warn(
                    f"Tag {dihedral_type.tag} not understood skipping"
                )
                continue
            children.append(Creator.parse_obj(dihedral_type.attrib))

        return cls(children=children)


class NonBondedAtom(ForceFieldChild):
    atom_type: str = Field(
        ..., alias="type", description="Atom Type Information"
    )

    charge: float = Field(
        ..., alias="charge", description="The charge of the atom type"
    )

    sigma: float = Field(..., alias="sigma", description="The Sigma parameter")

    epsilon: float = Field(
        ..., alias="epsilon", description="The epsilon parameter"
    )


@registers_loader(name="NonbondedForce")
class NonBondedForce(ForceFieldChild):
    children: List[NonBondedAtom]
    coulomb14scale: Optional[float] = None
    lj14scale: Optional[float] = None

    @classmethod
    def load_from_etree(cls, nonbonded_atoms):
        children = []
        coulomb14scale = nonbonded_atoms.attrib.get("coulomb14scale", None)
        lj14scale = nonbonded_atoms.attrib.get("lj14scale", None)

        for atom_type in nonbonded_atoms.iterchildren():
            if atom_type.tag == "Atom":
                children.append(NonBondedAtom.parse_obj(atom_type.attrib))

        return cls(
            children=children,
            coulomb14scale=coulomb14scale,
            lj14scale=lj14scale,
        )


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

    def iterate_on(self, children_type):
        if children_type not in loaders:
            raise ValueError(f"Only {list(loaders)} are supported")
        else:
            for child in self.children:
                if type(child) == loaders[children_type]:
                    for type_children in child.children:
                        yield type_children
