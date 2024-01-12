from setuptools import find_packages, setup

setup(
    name="forcefield_utilities",
    version="0.3.0",
    description="XML Conversion Utilities for MoSDeF ForceFields",
    authors="Umesh Timalsina",
    author_email="umesh.timalsina@vanderbilt.edu",
    install_requires=["foyer", "gmso"],
    license="MIT",
    packages=find_packages(),
    zip_safe=True,
    include_package_data=True,
)
