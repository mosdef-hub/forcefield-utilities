import os

from setuptools import find_packages, setup

setup(
    name="forcefield_utilities",
    version="0.0.1",
    description="XML Conversion Utilites for MoSDeF ForceFields",
    authors="Umesh Timalsina",
    author_email="umesh.timalsina@vanderbilt.edu",
    install_requires=["pydantic", "foyer", "lxml"],
    license="Apache-2.0",
    packages=find_packages(),
    zip_safe=True,
    include_package_data=True,
)
