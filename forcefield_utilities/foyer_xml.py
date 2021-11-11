import functools
import re
import warnings
from typing import Any, ClassVar, List, Optional, Set, Union

from gmso.core.angle_type import AngleType
from gmso.core.atom_type import AtomType
from gmso.core.bond_type import BondType
from gmso.core.dihedral_type import DihedralType
from gmso.core.forcefield import ForceField as GMSOForceField
from gmso.core.improper_type import ImproperType
from gmso.lib.potential_templates import (
    PotentialTemplate,
    PotentialTemplateLibrary,
)
from gmso.utils._constants import FF_TOKENS_SEPARATOR
from pydantic import BaseModel, Field

from forcefield_utilities.parameters_transformer import ParametersTransformer

__all__ = ["ForceField"]

GMSO_FF_WILDCARD = "*"

loaders = {}


class registers_loader:
    def __init__(self, name):
        self.name = name

    def __call__(self, child_class):
        functools.update_wrapper(self, child_class)
        loaders[self.name] = child_class
        return child_class


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

    def parameters(self):
        return self.dict(include={"length", "k"})


@registers_loader(name="HarmonicBondForce")
class HarmonicBondForce(ForceFieldChild):
    gmso_template: ClassVar[str] = "HarmonicBondPotential"

    children: List[Bond]

    def to_gmso_potentials(self, children):
        template = PotentialTemplateLibrary()[self.gmso_template]
        potentials = {"bond_types": {}}
        for child in self.children:
            parameters = ParametersTransformer.transform(
                self.gmso_template, child.parameters()
            )
            gmso_bond_type = BondType.from_template(template, parameters)

            if child.class1 and child.class2:
                gmso_bond_type.member_classes = [child.class1, child.class2]
                key = FF_TOKENS_SEPARATOR.join(gmso_bond_type.member_classes)
            else:
                gmso_bond_type.member_types = [child.type1, child.type2]
                key = FF_TOKENS_SEPARATOR.join(gmso_bond_type.member_types)

            potentials["bond_types"][key] = gmso_bond_type

        return potentials

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

    def parameters(self):
        return self.dict(include={"angle", "k"})


@registers_loader(name="HarmonicAngleForce")
class HarmonicAngleForce(ForceFieldChild):
    gmso_template: ClassVar[str] = "HarmonicAnglePotential"

    children: List[Angle] = []

    @classmethod
    def load_from_etree(cls, angles):
        children = []
        for angle_type in angles.iterchildren():
            if angle_type.tag == Angle.__name__:
                children.append(Angle.parse_obj(angle_type.attrib))
        return cls(children=children)

    def to_gmso_potentials(self, children):
        template = PotentialTemplateLibrary()[self.gmso_template]
        attr = "angle_types"
        potentials = {"angle_types": {}}
        for child in self.children:
            parameters = ParametersTransformer.transform(
                self.gmso_template, child.parameters()
            )
            gmso_angle_type = AngleType.from_template(template, parameters)

            if child.class1 and child.class2 and child.class3:
                gmso_angle_type.member_classes = [
                    child.class1,
                    child.class2,
                    child.class3,
                ]
                key = FF_TOKENS_SEPARATOR.join(gmso_angle_type.member_classes)
            else:
                gmso_angle_type.member_types = [
                    child.type1,
                    child.type2,
                    child.type3,
                ]
                key = FF_TOKENS_SEPARATOR.join(gmso_angle_type.member_types)

            potentials["angle_types"][key] = gmso_angle_type

        return potentials


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

    def parameters(self):
        return self.dict(include={"c0", "c1", "c2", "c3", "c4", "c5"})


class RBProper(RBDihedral):
    pass


class RBImproper(RBDihedral):
    pass


@registers_loader(name="RBTorsionForce")
class RBTorsionForce(ForceFieldChild):
    gmso_template: ClassVar[str] = "RyckaertBellemansTorsionPotential"

    children: List[Union[RBProper, RBImproper]]

    def to_gmso_potentials(self, children):
        template = PotentialTemplateLibrary()[self.gmso_template]
        potentials = {
            "dihedral_types": {},
            "improper_types": {},
        }
        count = 0
        for child in self.children:
            if isinstance(child, RBProper):
                Creator = DihedralType
                potential_dict = potentials["dihedral_types"]
            else:
                Creator = ImproperType
                potential_dict = potentials["improper_types"]

            parameters = ParametersTransformer.transform(
                self.gmso_template, child.parameters()
            )
            gmso_di_im_type = Creator.from_template(template, parameters)

            classes = [child.class1, child.class2, child.class3, child.class4]
            if any(classes):
                gmso_di_im_type.member_classes = [
                    class_ or "*" for class_ in classes
                ]
                key = FF_TOKENS_SEPARATOR.join(gmso_di_im_type.member_classes)
            else:
                gmso_di_im_type.member_types = [
                    child.type1,
                    child.type2,
                    child.type3,
                    child.type4,
                ]
                key = FF_TOKENS_SEPARATOR.join(gmso_di_im_type.member_types)

            potential_dict[key] = gmso_di_im_type

        return potentials

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

    def xml_dict(self):
        max_count = len(self.periodicity)
        xml_dict = {}

        for j in range(1, max_count + 1):
            xml_dict[f"periodicity{j}"] = self.periodicity[j - 1]
            xml_dict[f"phase{j}"] = self.phase[j - 1]
            xml_dict[f"k{j}"] = self.k[j - 1]

        return xml_dict

    def parameters(self):
        return self.dict(include={"periodicity", "phase", "k"})

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

        for key in attrib:
            if not (
                key.startswith("periodicity")
                or key.startswith("phase")
                or key.startswith("k")
            ):
                attrib_dict[key] = attrib[key]

        return attrib_dict

    @classmethod
    def parse_obj(cls, obj: Any):
        dict_attribs = cls.periodic_attribs_to_list(obj)
        return super().parse_obj(dict_attribs)


class PeriodicProper(PeriodicDihedral):
    gmso_template: ClassVar[str] = "PeriodicTorsionPotential"


class PeriodicImproper(PeriodicDihedral):
    gmso_template: ClassVar[str] = "PeriodicImproperPotential"


@registers_loader(name="PeriodicTorsionForce")
class PeriodicTorsionForce(ForceFieldChild):
    gmso_template: ClassVar[str] = "PeriodicTorsionPotential"

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

    def to_gmso_potentials(self, children):
        template = PotentialTemplateLibrary()[self.gmso_template]
        potentials = {
            "dihedral_types": {},
            "improper_types": {},
        }
        count = 0
        for child in self.children:
            if isinstance(child, PeriodicProper):
                Creator = DihedralType
                potential_dict = potentials["dihedral_types"]
            else:
                Creator = ImproperType
                potential_dict = potentials["improper_types"]

            parameters = ParametersTransformer.transform(
                self.gmso_template, child.parameters()
            )
            gmso_di_im_type = Creator.from_template(template, parameters)

            classes = [child.class1, child.class2, child.class3, child.class4]
            if any(classes):
                gmso_di_im_type.member_classes = [
                    class_ or "*" for class_ in classes
                ]
                key = FF_TOKENS_SEPARATOR.join(gmso_di_im_type.member_classes)
            else:
                gmso_di_im_type.member_types = [
                    child.type1,
                    child.type2,
                    child.type3,
                    child.type4,
                ]
                key = FF_TOKENS_SEPARATOR.join(gmso_di_im_type.member_types)

            potential_dict[key] = gmso_di_im_type

        return potentials


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

    def parameters(self):
        return self.dict(include={"charge", "sigma", "epsilon"})


@registers_loader(name="NonbondedForce")
class NonBondedForce(ForceFieldChild):
    gmso_template: ClassVar[str] = "LennardJonesPotential"
    children: List[NonBondedAtom]
    coulomb14scale: Optional[float] = None
    lj14scale: Optional[float] = None

    def to_gmso_potentials(self, children):
        template = PotentialTemplateLibrary()[self.gmso_template]
        nonbonded = {"atom_types": {}}
        foyer_atom_types = list(
            filter(lambda c: isinstance(c, AtomTypes), children)
        ).pop()
        atom_name_to_class_map = {}
        for type_ in foyer_atom_types.children:
            atom_name_to_class_map[type_.name] = type_

        for child in self.children:
            parameters = ParametersTransformer.transform(
                self.gmso_template, child.parameters()
            )
            gmso_atom_type = AtomType.from_template(template, parameters)
            # ToDO: Add Missing Properties
            foyer_atomtype = atom_name_to_class_map[child.atom_type]
            gmso_atom_type.name = foyer_atomtype.name
            gmso_atom_type.atomclass = foyer_atomtype.atom_class

            # doi, overrides, defintion, description
            gmso_atom_type.doi = foyer_atomtype.doi
            gmso_atom_type.overrides = (
                set(foyer_atomtype.overrides.strip().split(","))
                if foyer_atomtype.overrides
                else set()
            )
            gmso_atom_type.definition = foyer_atomtype.smarts_def
            gmso_atom_type.description = foyer_atomtype.desc

            # mass, charge
            gmso_atom_type.mass = foyer_atomtype.mass
            gmso_atom_type.charge = child.charge
            gmso_atom_type.tags = {"element": foyer_atomtype.element}

            nonbonded["atom_types"][child.atom_type] = gmso_atom_type

        return nonbonded

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

    def to_gmso_ff(self):
        ff = GMSOForceField()
        ff_potentials = {}
        for child in self.children:
            if hasattr(child, "gmso_template") and hasattr(
                child, "to_gmso_potentials"
            ):
                potentials = child.to_gmso_potentials(self.children)
                for attr in potentials:
                    if attr in ff_potentials:
                        ff_potentials[attr].update(potentials[attr])
                    else:
                        ff_potentials[attr] = potentials[attr]

        for attr in ff_potentials:
            setattr(ff, attr, ff_potentials[attr])
        try:
            nonbonded_force = list(
                filter(lambda c: isinstance(c, NonBondedForce), self.children)
            ).pop()
            ff.scaling_factors = {
                "electrostatics14Scale": nonbonded_force.coulomb14scale,
                "nonBonded14Scale": nonbonded_force.lj14scale,
            }
        except TypeError:
            warnings.warn("No nonbonded forces found")

        ff.name = self.name
        ff.version = self.version

        return ff
