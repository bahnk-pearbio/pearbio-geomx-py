"""
Analyse GeoMx samples.
"""

from setuptools import setup, find_packages

setup(
    name="pearbio_geomx_py",
    version="0.0.0",
    description="Tools to analyse GeoMx data.",
    author="Nourdine Bah",
    author_email="nourdine@pearbio.com",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["click"],
    entry_points={
        "console_scripts": [
            "extract_counts=pearbio_geomx_py.scripts.extract_counts:main",
            "gsea_patients=pearbio_geomx_py.scripts.gsea_patients:main",
        ]
    },
)
