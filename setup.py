import os
import sys
from shutil import rmtree

from setuptools import setup, Command

NAME = 'tiptop'
DESCRIPTION = '' #TODO
URL = 'https://github.com/FabioRossiArcetri/TIPTOP'
EMAIL = 'fabio.rossi@inaf.it'
AUTHOR = 'Fabio Rossi, INAF Arcetri Adaptive Optics group'
LICENSE = '' #TODO
KEYWORDS = 'Adaptive Optics, Astrophysics, INAF, Arcetri',

here = os.path.abspath(os.path.dirname(__file__))
# Load the package's __version__.py module as a dictionary.
about = {}
with open(os.path.join(here, NAME, '__version__.py')) as f:
    exec(f.read(), about)

setup(name=NAME,
      description=DESCRIPTION,
      version=about['__version__'],
      classifiers=['Development Status :: 4 - Beta',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python :: 3',
                   ],
      long_description=open('README.md').read(),
      url=URL,
      author_email=EMAIL,
      author=AUTHOR,
      license=LICENSE,
      keywords=KEYWORDS,
      packages=['tiptop',
                ],
      install_requires=["numpy",
                        "scipy",
                        "matplotlib",
                        "astropy",
                        ],
      )
