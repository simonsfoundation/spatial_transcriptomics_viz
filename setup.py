
# Thanks: http://pythonhosted.org/an_example_pypi_project/setuptools.html

import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "spatial_transcriptomics_viz",
    version = "0.0.1",
    author = "Aaron Watters",
    author_email = "awatters@simonsfoundation.org",
    description = ("visualization aids for spatial transcriptomics."),
    license = "BSD",
    keywords = "jupyter widget genomics",
    url = "http://packages.python.org/spatial_transcriptomics_viz",
    packages=['spatial_transcriptomics_viz'],
    #package_data={'jp_gene_viz': ['*.js']},
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: BSD License",
    ],
)
