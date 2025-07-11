"""
Classes Module.

This module contains all the main classes for the transmission_models package.

Main Classes
------------
host : Host class representing infected individuals
didelot_unsampled : Main class implementing the Didelot et al. (2017) framework
genetic_prior_tree : Prior distribution for genetic sequence data
location_distance_prior_tree : Prior distribution for location distance data
same_location_prior_tree : Prior distribution for same location probability
MCMC : Markov Chain Monte Carlo sampling algorithms

Submodules
----------
mcmc : MCMC sampling classes and algorithms
"""

from .host import host
from .didelot_unsampled import didelot_unsampled
from .genetic_prior import genetic_prior_tree
from .location_prior import location_distance_prior_tree, same_location_prior_tree
from .mcmc import MCMC

__all__ = [
    'host',
    'didelot_unsampled', 
    'genetic_prior_tree',
    'location_distance_prior_tree',
    'same_location_prior_tree',
    'MCMC'
]
