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

The Jupyter Notebook [Example.ipynb](https://github.com/oscarf/transmission_models/blob/main/Example.ipynb) contains an illustration of the usage of this library. 

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

See the [Example.ipynb](https://github.com/oscarcapote/transmission_models/blob/main/examples/Example.ipynb) notebook for detailed usage examples and the [documentation](https://www.maths.usyd.edu.au/u/oscarf/transmission_models_documentation/) for complete API reference.

# Tree Visualization JS library:

`tree_plot.js` is an interactive tree visualization library build in JavaScript and D3.js. This tool allows you to:

- **Interactive visualization** of transmission trees with D3.js
- **Customizable node colors** for sampled and unsampled hosts
- **Tooltips** showing host attributes on hover
- **Toggle between layouts**: classic tree layout and infection time-based layout
- **Responsive design** that adapts to window resizing

## Example Visualization

In https://www.maths.usyd.edu.au/u/oscarf/tree_layout/ you can upload your jsons to visualize your sampled networks:

![Tree Layout Visualization](examples/Screenshot%202025-07-30%20at%2015-49-21%20Tree%20layout.png)

This webpage have been developed using tree_layout.js

*Example of a tree layout visualization generated using the interactive tree visualization tool.*

## Documentation of tree_layout.js

For complete documentation on using the tree visualization features, see the [Tree Visualization Guide](./docs/tree_visualization.md) in the documentation.

## Online Tool

You can also use the online tree layout tool at: https://www.maths.usyd.edu.au/u/oscarf/tree_layout/
