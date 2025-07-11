"""
Transmission Models Package.

This package provides tools for modeling viral transmission networks
using phylogenetic and epidemiological data. It implements the Didelot
et al. (2017) framework for transmission tree inference with unsampled
hosts by using a MCMC sampling.

Main modules:
    - classes: Core classes including the Didelot unsampled model, host class, and priors
    - utils: Utility functions for tree manipulation and visualization

The package supports:
    - Bayesian inference using MCMC sampling
    - Integration of genetic sequence data
    - Location-based transmission modeling
    - Visualization of transmission networks

References
----------
Didelot, X., Gardy, J., & Colijn, C. (2017). Bayesian inference of
transmission chains using timing of events, contact and genetic data.
PLoS computational biology, 13(4), e1005496.
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from random import choice,randint,random,sample,choices
from scipy.stats import nbinom, gamma, binom, expon, norm
from matplotlib.lines import Line2D
import pandas as pd
from networkx.drawing.nx_pydot import graphviz_layout

# Importing classes
from .classes import host, didelot_unsampled, genetic_prior_tree, location_distance_prior_tree, same_location_prior_tree, MCMC



# Importing utilities
from . import utils

# Importing necessary libraries
import scipy.special as sc
import scipy.stats as st
import numpy as np


__all__ = [
    'host',
    'didelot_unsampled',
    'genetic_prior_tree',
    'location_distance_prior_tree',
    'same_location_prior_tree',
    'MCMC',
]