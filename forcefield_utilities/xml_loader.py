import importlib
import os
import pathlib

from lxml import etree

from forcefield_utilities.foyer_xml import ForceField


def call_on_import(func):
    """Declare a decorator that will run `func` when imported."""
    func()


def _get_xml_path(from_package, relative_path):
    """Use source of a python package to locate and cache the address of a file."""
    from pkg_resources import resource_filename

    return resource_filename(from_package, relative_path)


def _get_foyer_forcefield(xml_path):
    """Return the foyer Forcefield object from the relative path ``xml_path`` inside the foyer package."""
    with open(xml_path) as ff_file:
        root = etree.parse(ff_file).getroot()
        return ForceField.load_from_etree(root)


def load_xml(xml_name, rel_to_module=False):
    """Return the foyer Forcefield object for the `xml_name file found in working directory or in foyer.

    Parameters
    __________
    xml_name : str
        Additional name of xml files in the foyer library
    rel_to_module : bool, optional default=False
        A flag to determine where to start lookup for path. If false,
        start from current working directory. If true, lookup inside
        foyer package at foyer/forcefields/xml/
    """
    try:
        if not rel_to_module:
            ff = _get_foyer_forcefield(
                _get_xml_path("foyer", f"forcefields/xml/{xml_name}")
            )
        else:
            ff = _get_foyer_forcefield(xml_name)
    except FileNotFoundError as e:
        if not os.path.splitext(xml_name)[1] == ".xml":
            print_error_message = (
                r"Please make sure the path to the xml file, "
                "uses the proper XML extension (.xml)."
            )
            raise ValueError(print_error_message)
        else:
            raise e
    return ff


def load_oplsaa():
    """Return the foyer Forcefield object for the oplsaa.xml file found in the foyer/forcefields/xml/ directory."""
    return _get_foyer_forcefield(
        _get_xml_path("foyer", "forcefields/xml/oplsaa.xml")
    )


def load_trappe_ua():
    """Return the foyer Forcefield object for the trappe-ua.xml file found in the foyer/forcefields/xml/ directory."""
    return _get_foyer_forcefield(
        _get_xml_path("foyer", "forcefields/xml/trappe-ua.xml")
    )


class FoyerFFs:
    """Object to provide methods to forcefields shipped with Foyer.

    Attributes
    __________
    ff_registry : dict, keys are strings and values are methods to load a given xml
        This is a default dictionary which directs the method to given default loaders.
    loaded_ffs : dict, keys are strings, values are the loaded forcefield
        This is a place to store loaded xmls, which can be accessed through
        the getter by indexing the FoyerFFs object. xmls are only stored once,
        and custom_xml is where a forcefield is stored from path.

    Methods
    __________
    get_ff(ffname="oplsaa"):
        Load and directly return the forcefield xml from Foyer or a path.
        Parameters
        __________
        ffname : str
            Name of forcefield to load. See self.ff_registry attribute
            for loading method.
    load(ffname="oplsaa", rel_to_module=False):
        Load and return the forcefield xml from Foyer or a path.
    """

    ff_registry = {
        "oplsaa": load_oplsaa,
        "trappe_ua": load_trappe_ua,
        "custom_xml": load_xml,
    }

    def __init__(self):
        super().__init__()
        self.loaded_ffs = {}

    @classmethod
    def get_ff(self, ffname, rel_to_module=False):
        """Load and directly return the forcefield xml from Foyer or a path.
        Parameters
        __________
        ffname : str
            Name of forcefield to load. See self.ff_registry attribute
            for loading method.
        rel_to_module : bool, optional default=False
            A flag to determine where to start lookup for path. If false,
            start from current working directory. If true, lookup inside
            foyer package at foyer/forcefields/xml/

        Returns
        ______
        ff : forcefield_utilities.foyer_xml.Forcefield
        """
        if ffname in self.ff_registry.keys():
            ff = self.ff_registry[ffname]()
        else:
            ff = self.ff_registry["custom_xml"](ffname, rel_to_module)
        return ff

    def load(self, ffname, rel_to_module=False):
        """Load and return the forcefield xml from Foyer or a path.
        Parameters
        __________
        ffname : str
            Name of forcefield to load. See self.ff_registry attribute
            for loading method.
        rel_to_module : bool, optional default=False
            A flag to determine where to start lookup for path. If false,
            start from current working directory. If true, lookup inside
            foyer package at foyer/forcefields/xml/

        Returns
        ______
        ff : forcefield_utilities.foyer_xml.Forcefield
        """
        if ffname in ["oplsaa", "trappe_ua"]:
            self.loaded_ffs[ffname] = self.ff_registry[ffname]()
        else:
            self.loaded_ffs[ffname] = self.ff_registry["custom_xml"](
                ffname, rel_to_module
            )
        return self.loaded_ffs[ffname]

    def __getitem__(self, ffname):
        """Get function for indexing by loaded forcefields"""
        if ffname not in self.loaded_ffs.keys():
            return None
        else:
            return self.ff_registry[ffname]()


@call_on_import
def register_gaff():
    """Include GAFF as part of FoyerFFs if antefoyer is installed locally."""
    try:
        importlib.import_module("antefoyer")
    except ImportError:
        return
    gaff_path = _get_xml_path("antefoyer", "xml/gaff.xml")
    FoyerFFs.gaff = _get_foyer_forcefield(gaff_path)
