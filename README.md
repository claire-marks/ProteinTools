# ProteinTools - README

### Overview

ProteinTools is a collection of useful (to me, anyway!) scripts and functions for dealing with proteins. It is written in Python 3 and consists of a module and a number of scripts for more specific jobs.

Who am I? I'm a Research Software Engineer at the University of Oxford, working in the Oxford Protein Informatics Group (OPIG). I've done a lot of research in protein loop modelling, hence the large number of loop-related code found here!

This is by no means a complete piece of work; essentially it exists to make my life easier during day-to-day research and so I will continue to add to it!

### Requirements

The following Python packages are required:
  - numpy
  - Biopython
  - requests

Also required:
  - [DSSP](https://github.com/cmbi/dssp) - the DSSP executable should be in your PATH.


### Installation

Installation should be as easy as this:

`python setup.py install`

Or:

`python setup.py install --user`


### The 'protmod' Module

I've tried to group similarly-themed functions together:

  - protmod = basic functions such as converting a one-letter amino acid code to a three-letter one
  - protmod.webtools = functions for interacting with the web, e.g. downloading files from the PDB
  - protmod.structures = functions for dealing with protein structures, normally using the Biopython PDB parser - for example, calculating phi/psi angles, selecting fragments
  - protmod.loops = specific functions involving loops
  - protmod.sequences = 


### The Scripts


