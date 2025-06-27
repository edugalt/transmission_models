# Transmission Model Library

This repository contains the codes used in the manuscript "Inference of epidemic networks: the effect of different data
types", by Oscar Fajardo-Fontiveros, Carl J. E. Suster, and Eduardo G. Altmann, developed at The University
of Sydney, Sydney, NSW, Australia

# Repository description

The Jupyter Notebook "Example.ipynb" contains an illustration of the usage of this library. 

- data/

Contains the data used in the manuscript (49 cases of positive tests in NSW, from mid 2021). Information about the time, genetic distance, and location of each case is provided in separate files.

- models/

Contains the code that implements the transmission model and the MCMC used to sample transmission trees.

- priors/

Contains the code that implements the genetic and location models.

- utils/

Contains additional functions used in the library.
