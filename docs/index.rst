Transmission Models
===================

Python Version Documentation

A Python library for modeling viral transmission networks using phylogenetic and epidemiological data. This package implements the Didelot et al. (2017) framework for transmission tree inference with unsampled hosts.

.. toctree::
   :maxdepth: 5
   :caption: Contents:

   classes
   functions
   tree_visualization


Classes
=======

See the :doc:`classes` page for detailed class documentation.

Functions
=========

See the :doc:`functions` page for a complete list of functions and their usage.

Tree Visualization
==================

Interactive visualization of transmission trees is supported via a JavaScript library. See the :doc:`tree_visualization` for full documentation, usage examples, and API reference.

Overview
--------

This library provides tools for Bayesian inference of transmission chains using timing of events, contact and genetic data. It implements the Mixed-Membership Stochastic Block Model approach to transmission network analysis, allowing researchers to:

* Reconstruct transmission networks from genetic and epidemiological data
* Infer unsampled hosts in transmission chains
* Incorporate genetic sequence evolution in transmission modeling
* Perform Bayesian inference using MCMC sampling
* Visualize transmission networks and phylogenetic relationships

Installation
-----------

Prerequisites
~~~~~~~~~~~~

This package requires Python 3.8 or higher and the following dependencies:

* numpy
* pandas
* scipy
* matplotlib
* seaborn

Installing the Package
~~~~~~~~~~~~~~~~~~~~~

To install the package in editable mode for development:

.. code-block:: bash

   git clone https://github.com/oscarcapote/transmission_models.git
   cd transmission_models
   pip install -e .

This will install the package in editable mode, meaning any changes you make to the source code will be immediately available without reinstalling.

Verifying Installation
^^^^^^^^^^^^^^^^^^^^^

After installation, you can verify the package is working correctly:

.. code-block:: python

   from transmission_models import *
   print("Transmission Models package installed successfully!")

Features
--------

* **Transmission Network Modeling**: Build and analyze transmission networks with unsampled hosts
* **Phylogenetic Integration**: Incorporate genetic sequence data.
* **MCMC Sampling**: Bayesian inference using Markov Chain Monte Carlo methods.
* **Prior Distributions**: Flexible prior specification for genetic and location data
* **Visualization**: Tools for plotting transmission networks and phylogenetic trees
* **High-performance**: Optimized algorithms for large-scale transmission analysis


Requirements
------------

* Python >= 3.8
* Required packages:
  * numpy
  * scipy
  * networkx
  * matplotlib
  * pandas
* Optional but recommended:
  * imageio (for animation support)

Usage
-----

The library can be used for transmission network inference where you have:
* Genetic sequence data from sampled hosts
* Epidemiological data (sampling times, infection times)
* Optional location or contact data
* Need to infer unsampled hosts in transmission chains

Basic Example
-------------

.. code-block:: python

   from transmission_models import host
   from transmission_models.models import didelot_unsampled
   
   # Define model parameters
   sampling_params = {
       "pi": 0.1,        # sampling probability
       "k_samp": 2.0,    # shape parameter for gamma distribution
       "theta_samp": 1.0 # scale parameter for gamma distribution
   }
   
   offspring_params = {
       "r": 1.5,         # rate of infection
       "p_inf": 0.3      # probability of infection
   }
   
   infection_params = {
       "k_inf": 2.0,     # shape parameter for gamma distribution
       "theta_inf": 1.0  # scale parameter for gamma distribution
   }
   
   # Create a transmission model
   model = didelot_unsampled(sampling_params, offspring_params, infection_params)
   
   # Add a root host
   root = model.add_root(t_sampl=10, id="root", genetic_data=['A', 'T', 'C', 'G'], t_inf=0)
   
   # Run MCMC sampling
   from transmission_models.models.MCMC import MCMC
   mcmc = MCMC(model)
   
   for i in range(1000):
       move, gg, pp, P, accepted, DL = mcmc.MCMC_iteration()
       if i % 100 == 0:
           print(f"Iteration {i} - Move: {move}, Accepted: {accepted}")

How It Works
------------

The Transmission Models library implements the Didelot et al. (2017) framework for transmission tree inference. This approach:

* Models transmission networks as directed trees
* Incorporates genetic sequence evolution using mutation models
* Handles unsampled hosts through Bayesian inference
* Uses MCMC sampling to explore the posterior distribution

The model combines three main components:

1. **Sampling Model**: Gamma distribution for sampling times
2. **Offspring Model**: Negative binomial distribution for offspring number
3. **Infection Model**: Gamma distribution for infection times

Documentation and Examples
--------------------------

* Download the package from the repository
* Tutorial: See Example.ipynb for detailed usage examples
* API Reference: Complete documentation of all classes and functions

References
----------

1. Didelot, X., Gardy, J., & Colijn, C. (2017). Bayesian inference of transmission chains using timing of events, contact and genetic data. PLoS computational biology, 13(4), e1005496.

License
-------

This project is licensed under the MIT License - see the LICENSE file for details.

Author
------

Transmission Models Team

Indices and tables
===================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search` 