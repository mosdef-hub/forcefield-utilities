from lxml import etree

from forcefield_utilities.foyer_xml import ForceField


def get_foyer_forcefield(filename):
    from pkg_resources import resource_filename

    xml_path = resource_filename("foyer", f"forcefields/xml/{filename}")
    with open(xml_path) as ff_file:
        root = etree.parse(ff_file).getroot()
        return ForceField.load_from_etree(root)


class FoyerFFs:
    oplsaa = get_foyer_forcefield("oplsaa.xml")
