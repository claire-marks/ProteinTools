#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

numpy_include_dir = np.get_include()

setup(name = "protmod",
      version = "1.0",
      description = "Useful things for dealing with proteins!",
      author = "Claire Marks",
      author_email = "marks@stats.ox.ac.uk",
      packages = [
          "protmod",
      ],
      package_dir = {
          "protmod": "lib/python/protmod",
      },
      scripts = ["bin/LoopFinder", "bin/DistanceFromGermline"]
)
