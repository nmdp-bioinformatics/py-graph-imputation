'''
Created on Feb 16, 2018

@author: mmaiers-nmdp
'''

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='imputegl',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.0.4',

    description='Imputation based on Neo4j graph',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/nmdp-bioinformatics/graph-imputation-match',

    packages=[
        'imputegl'
    ],
    package_dir={'imputegl':
                 'imputegl'},
    # Author details
    author_email='mmaiers@nmdp.org',

    # Choose your license
    license='LGPL3',

    # What does your project relate to?
    keywords='imputation HLA graph',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['pandas>=1.0.3', 'networkx==2.3'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },

)


