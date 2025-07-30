Transmission Models
===================

Python Version Documentation

A Python library for modeling viral transmission networks using phylogenetic and epidemiological data. This package implements the Didelot et al. (2017) framework for transmission tree inference with unsampled hosts.

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

* Download the package from the `repository <https://github.com/oscarcapote/transmission_models>`_
* Tutorial: See :doc:`tutorial` for a complete workflow example with plots
* Jupyter Notebook: See `Example.ipynb <https://github.com/oscarcapote/transmission_models/blob/main/examples/Example.ipynb>`_ for the original notebook
* API Reference: Complete documentation of all classes and functions

References
----------

1. Didelot, X., Gardy, J., & Colijn, C. (2017). Bayesian inference of transmission chains using timing of events, contact and genetic data. PLoS computational biology, 13(4), e1005496.



Indices and tables
===================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search` 

.. toctree::
   :maxdepth: 5
   :caption: Contents:

   tutorial
   classes
   functions
   tree_visualization


