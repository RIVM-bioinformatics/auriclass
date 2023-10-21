import sys

from setuptools import find_packages, setup

from auriclass.version import __version__

if sys.version_info.major != 3:
    print("Error: you must execute setup.py using Python 3")
    sys.exit(1)

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="AuriClass",
    description="Fast prediction of Candida auris WGS data",
    author="Boas van der Putten",
    author_email="ids-bioinformatics@rivm.nl",
    license="AGPLv3",
    version=__version__,
    packages=find_packages(),
    python_requires=">=3",
    install_requires=[
        "pandas",
        "pyfastx",
    ],
    entry_points={
        "console_scripts": [
            "auriclass = auriclass.main:main",
        ]
    },
    keywords=[],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
