# Transmission Model Library

This repository contains the codes and data used in the manuscript "Inference of epidemic networks: the effect of different data
types", by [Oscar Fajardo-Fontiveros](https://www.maths.usyd.edu.au/u/oscarf/), [Carl J. E. Suster](https://www.sydney.edu.au/medicine-health/about/our-people/academic-staff/carl.suster.html), and [Eduardo G. Altmann](https://www.maths.usyd.edu.au/u/ega), developed at The University
of Sydney, Sydney, NSW, Australia. The code in this repository was created by [Oscar Fajardo-Fontiveros](https://www.maths.usyd.edu.au/u/oscarf/).

# Installation

## Prerequisites

This package requires Python 3.8 or higher and the following dependencies:
- numpy
- pandas
- scipy
- matplotlib
- seaborn

## Installing the Package

To install the package in editable mode for development:

```bash
git clone https://github.com/oscarcapote/transmission_models.git
cd transmission_models
pip install -e .
```

This will install the package in editable mode, meaning any changes you make to the source code will be immediately available without reinstalling.

## Verifying Installation

After installation, you can verify the package is working correctly:

```python
from transmission_models import *
print("Transmission Models package installed successfully!")
```

# Repository description

The Jupyter Notebook [Example.ipynb](https://github.com/edugalt/transmission_models/blob/main/Example.ipynb) contains an illustration of the usage of this library. 

- data/

Contains the data used in the manuscript (49 cases of positive tests in NSW, from mid 2021). Information about the time, genetic distance, and location of each case is provided in separate files.

- src/transmission_models/

Contains the main package code organized as follows:

- classes/

Contains the code that implements the transmission model and the MCMC used to sample transmission trees.

- priors/

Contains the code that implements the genetic and location models.

- utils/

Contains additional functions used in the library.

# Usage

See the [Example.ipynb](https://github.com/edugalt/transmission_models/blob/main/Example.ipynb) notebook for detailed usage examples and the [documentation](https://oscarcapote.github.io/transmission_models/) for complete API reference.
