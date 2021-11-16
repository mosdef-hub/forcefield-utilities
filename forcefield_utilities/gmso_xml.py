from typing import List, Optional, Union

from pydantic import BaseModel, Field


class GMSOXMLTag(BaseModel):
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

    value: float = Field(
        ..., description="The value of the parameter", alias="value"
    )

    @classmethod
    def load_from_etree(cls, root):
        return cls(**root.attrib)


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
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "Parameters":
                children.append(Parameters.load_from_etree(el))
        return cls(children=children, **root.attrib)


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

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "ParametersUnitDef":
                children.append(ParametersUnitDef.load_from_etree(el))
            elif el.tag == "BondType":
                children.append(BondType.load_from_etree(el))
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
        attribs = root.attrib
        children = []
        for el in root.iterchildren():
            if el.tag == "Parameters":
                children.append(Parameters.load_from_etree(el))
        return cls(children=children, **root.attrib)


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
        description="Type 3 for this Dihedral/Improper type",
        alias="type3",
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
        return cls(children=children, **root.attrib)


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
            children.append(Parameter.load_from_etree(el))
        return cls(children=children, **attribs)


class PairPotentialTypes:
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

    @staticmethod
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

    electrostatic_14_scale: float = Field(0.5, alias="electrostatics14Scale")

    nonbonded_14_scale: float = Field(0.5, alias="nonBonded14Scale")

    @classmethod
    def load_from_etree(cls, root):
        attribs = root.attrib
        children = []
        for unit in root.iterchildren():
            children.append(Units.load_from_etree(unit))
        return cls(children=children, **attribs)


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
        return cls(children=children, **attribs)
