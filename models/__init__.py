"""
Models Module.

This module contains the core transmission modeling classes and algorithms
for viral transmission network inference.

Main Classes
------------
didelot_unsampled : Main class implementing the Didelot et al. (2017) framework
                   for transmission tree inference with unsampled hosts.

Submodules
----------
MCMC : Markov Chain Monte Carlo sampling algorithms
topology_movements : Functions for modifying transmission tree topology

References
----------
Didelot, X., Gardy, J., & Colijn, C. (2017). Bayesian inference of
transmission chains using timing of events, contact and genetic data.
PLoS computational biology, 13(4), e1005496.
"""

from transmission_models.models import *
from transmission_models.models.didelot_unsampled import *
from transmission_models.models.MCMC import *
from transmission_models.utils import *
# from transmission_models.utils import tree_to_newick
