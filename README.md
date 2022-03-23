[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/mosdef-hub/forcefield-utilities/main.svg)](https://results.pre-commit.ci/latest/github/mosdef-hub/forcefield-utilities/main)
[![test](https://github.com/mosdef-hub/forcefield-utilities/actions/workflows/CI.yaml/badge.svg)](https://github.com/mosdef-hub/forcefield-utilities/actions/workflows/CI.yaml)

## forcefield-utilities
This repository contains utility functions and classes necessary to convert xml forcefields to `foyer` and `gmso` Forcefields.

### Install Instructions
This package is available on conda-forge. Use the following command to install:
```bash
$ conda install -c conda-forge forcefield-utilities
```

How to use this package to load foyer or GMSO forcefields:
```python
import forcefield_utilities as ff_utils

#get foyer and gmso forcefields stored in foyer/forcefields
foyer_opls = ff_utils.FoyerFFs.GetFF('oplsaa')
gmso_opls = foyer_opls.to_gmso_ff()

#get foyer and gmso forcefields from a local xml
ffxml_loader = ff_utils.FoyerFFs()
my_foyer_xml = ffxml_loader.load("PATH_TO_FILE/MYFILE.xml", rel_to_module=False)
my_gmso_forcefield = my_foyer_xml['MYFILE.xml'].to_gmso_ff()
```
