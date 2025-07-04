"""
Transmission Models Package.

This package provides tools for modeling viral transmission networks
using phylogenetic and epidemiological data. It implements the Didelot
et al. (2017) framework for transmission tree inference with unsampled
hosts by using a MCMC sampling.

Main modules:
    - models: Core transmission modeling classes including the Didelot unsampled model
    - priors: Prior probability distributions for genetic and location data
    - utils: Utility functions for tree manipulation and visualization
    - host: Host class representing infected individuals

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




# Importing models
from transmission_models.models.didelot_unsampled import didelot_unsampled
from transmission_models.host import host

# Importing priors
from transmission_models.priors.genetic_prior import *
from transmission_models.priors.location_prior import *

# Importing utilities
import transmission_models.utils as utils

# Importing necessary libraries
import scipy.special as sc
import scipy.stats as st
import numpy as np
