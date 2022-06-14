import importlib
import os
from pathlib import Path
from typing import Union

from lxml import etree

from forcefield_utilities.foyer_xml import ForceField
from forcefield_utilities.utils import (
    call_on_import,
    deprecate_kwargs,
    get_package_file_path,
)

custom_forcefields = {}


def _parse_foyer_xml(xml_path):
    """Return the foyer Forcefield object from the relative path ``xml_path`` inside the foyer package."""
    with open(xml_path) as ff_file:
        root = etree.parse(ff_file).getroot()
        return ForceField.load_from_etree(root)


class FoyerFFs:
    """Object to provide methods to forcefields shipped with Foyer.

    Attributes
    __________
    loaded_ffs : dict, keys are strings, values are the loaded forcefield
        This is a place to store loaded xmls, which can be accessed through
        the getter by indexing the FoyerFFs object. xmls are only stored once,
        and custom_xml is where a forcefield is stored from path. Note: This is a
        class level attribute.

    Methods
    _______
    get_ff(ffname="oplsaa"):
        Load and directly return the forcefield xml from Foyer or a path.
        Parameters
        __________
        ffname : str
            Name of forcefield to load. See self.ff_registry attribute
            for loading method.
    load(ffname="oplsaa", rel_to_module=False):
        Load and return the forcefield xml from Foyer or a path.

    Notes
    -----
    This class caches it caches the loaded forcefield files. A method called
    `clear_loaded_ffs` can be used to clear the cache.
    """

    loaded_ffs = {}
    overwritten_custom_ffs = set()

    @classmethod
    @deprecate_kwargs(deprecated_kwargs={"rel_to_module"})
    def get_ff(
        cls, ffname: Union[str, Path], rel_to_module: bool = False
    ) -> ForceField:
        """Load and directly return the forcefield xml from Foyer or a path.
        Parameters
        __________
        ffname : str, pathlib.Path
            Name of forcefield to load or path of the foyer XML file
        rel_to_module : bool, optional default=False
            Deprecated, has no effect

        Returns
        _______
        ff : forcefield_utilities.foyer_xml.Forcefield
        """
        loader = cls()
        return loader.load(ffname, rel_to_module)

    @deprecate_kwargs(deprecated_kwargs={"rel_to_module"})
    def load(self, ffname, rel_to_module=False):
        """Load and return the forcefield xml from Foyer or a path.
        Parameters
        __________
        ffname : str, pathlib.Path
              Name of forcefield to load or path of the foyer XML file
        rel_to_module : bool, optional default=False
             Deprecated, has no effect

        Returns
        ______
        ff : forcefield_utilities.foyer_xml.Forcefield

        Notes
        -----
        This method returns a parsed foyer XML by either taking the path of the foyer XML or If ffname is a string that doesn't end in `.xml`,
        it will check the `custom_forcefields` global dictionary for the path and if not found
        check the foyer/forcefields/xml directory for a file named `ffname`+'.xml' and try to
        load it from there.

        Warnings
        --------
        This method caches the loaded forcefields by name, please use the clear_cache method
        if you changed the forcefield XML file. This might be especially relavent when using from a jupyter notebook.
        """
        if (
            ffname in self.loaded_ffs
            and ffname not in self.overwritten_custom_ffs
        ):
            return self.loaded_ffs[ffname]

        ff_path = Path(ffname)

        if ffname in custom_forcefields:
            self.loaded_ffs[ffname] = _parse_foyer_xml(
                xml_path=custom_forcefields[ffname]
            )
            self.overwritten_custom_ffs.discard(ffname)
            return self.loaded_ffs[ff_path.name]

        if self._is_xml(ff_path):
            self.loaded_ffs[ff_path.name] = _parse_foyer_xml(
                xml_path=ff_path.resolve()
            )
        else:
            xml_path = get_package_file_path(
                "foyer", f"forcefields/xml/{ffname}.xml"
            )
            self.loaded_ffs[ffname] = _parse_foyer_xml(xml_path)

        return self.loaded_ffs[ffname]

    def __getitem__(self, ffname) -> ForceField:
        """Get function for indexing by loaded forcefields."""
        return self.loaded_ffs[ffname]

    def __getattr__(self, item) -> ForceField:
        """Accessor for loaded forcefields."""
        if item in self.loaded_ffs:
            return self.loaded_ffs[item]
        else:
            raise AttributeError(
                f"{self.__class__.__name__} has no attribute {item}"
            )

    def register_custom_forcefield(
        self, name: str, path_: Union[str, Path], overwrite: bool = True
    ) -> None:
        """Register a custom foyer forcefield's XML path to load.

        Parameters
        ----------
        name: str
            The unique name of the custom forcefield
        path_: str, Path
            The XML path for the forcefield XML file
        """
        if not overwrite and name in custom_forcefields:
            raise ValueError(
                f"Forcefield {name} is already registered to point to "
                f"{custom_forcefields[name]}. Please use overwrite=True "
                f"if you wish to overwrite it."
            )
        if overwrite and name in custom_forcefields:
            self.overwritten_custom_ffs.add(name)

        custom_forcefields[name] = str(path_)

    def clear_cache(self):
        self.loaded_ffs = {}
        self.overwritten_custom_ffs = set()

    @staticmethod
    def _is_xml(path_: Path) -> bool:
        return path_.suffix == ".xml"


@call_on_import
def register_gaff():
    """Include GAFF as part of FoyerFFs if antefoyer is installed locally."""
    try:
        importlib.import_module("antefoyer")
    except ImportError:
        return
    custom_forcefields["gaff"] = get_package_file_path(
        "antefoyer", "xml/gaff.xml"
    )
