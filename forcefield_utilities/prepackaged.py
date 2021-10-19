import importlib

from lxml import etree

from forcefield_utilities.foyer_xml import ForceField


def call_on_import(func):
    func()


def get_xml_path(from_package, relative_path):
    from pkg_resources import resource_filename

    return resource_filename(from_package, relative_path)


def get_foyer_forcefield(xml_path):
    with open(xml_path) as ff_file:
        root = etree.parse(ff_file).getroot()
        return ForceField.load_from_etree(root)


def load_oplsaa():
    return get_foyer_forcefield(
        get_xml_path("foyer", "forcefields/xml/oplsaa.xml")
    )


def load_trappe_ua():
    return get_foyer_forcefield(
        get_xml_path("foyer", "forcefields/xml/trappe-ua.xml")
    )


class FoyerFFs:
    oplsaa = load_oplsaa()
    trappe_ua = load_trappe_ua()


@call_on_import
def register_gaff():
    try:
        importlib.import_module("antefoyer")
    except ImportError:
        return
    gaff_path = get_xml_path("antefoyer", "xml/gaff.xml")
    FoyerFFs.gaff = get_foyer_forcefield(gaff_path)
