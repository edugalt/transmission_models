Tree Visualization Library Documentation
=======================================

Overview
--------

``tree_plot.js`` is a JavaScript library for visualizing transmission tree structures using D3.js. It supports interactive features like tooltips to see host attributes and color customization.

Installation
------------

1. **Include D3.js** in your HTML::

   <script src="https://d3js.org/d3.v6.min.js"></script>

2. **Include ``tree_plot.js``** in your HTML::

   <script src="https://www.maths.usyd.edu.au/u/oscarf/tree_layout/tree_plot.js"></script>

Usage
-----

1. Prepare your HTML
~~~~~~~~~~~~~~~~~~~~
Add a container ``<div id="myTreeContainer"></div>`` where the tree will be rendered::

   <div id="myTreeContainer"></div>

2. Prepare your data
~~~~~~~~~~~~~~~~~~~~
Your data should be in JSON format, representing a tree structure. The library supports two formats:

**Format 1: Direct tree object**::

   {
     "name": "Root",
     "Infection time": 0,
     "Sampled": true,
     "Sampling time": 5,
     "children": [
       {
         "name": "Child 1",
         "Infection time": 1,
         "Sampled": false,
         "children": []
       }
     ]
   }

**Format 2: Tree object inside a 'tree' attribute**::

   {
     "tree": {
       "name": "Root",
       "Infection time": 0,
       "Sampled": true,
       "Sampling time": 5,
       "children": [
         {
           "name": "Child 1",
           "Infection time": 1,
           "Sampled": false,
           "children": []
         }
       ]
     }
   }

The library will automatically detect and use the correct format.

3. Call ``createTreeMap``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: javascript

   createTreeMap("treeData.json", "#myTreeContainer", {
     sampledColor: "#1f77b4", // Optional: color for sampled nodes
     unsampledColor: "#ff7f0e", // Optional: color for unsampled nodes
     tooltip: customTooltip, // Optional: custom tooltip function or d3 selection
     toggleSwitchId: "mySwitch" // Optional: custom toggle switch ID
   });

API Reference
-------------

``createTreeMap(jsonFilename, containerDiv, options)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- **``jsonFilename``** (string): Path to your tree data JSON file.
- **``containerDiv``** (string): CSS selector for the container div (e.g., ``#myTreeContainer``).
- **``options``** (object, optional):

  - **``sampledColor``** (string): Color for sampled nodes (default: ``"dodgerblue"``).
  - **``unsampledColor``** (string): Color for unsampled nodes (default: ``"darkorange"``).
  - **``tooltip``** (function or d3 selection): Custom tooltip handler.
    - If a function, it will be called as ``tooltip(event, d, html_content, isMouseOut)``.
    - If a d3 selection, it will be used as the tooltip div.
  - **``toggleSwitchId``** (string): ID of the toggle switch element (default: ``"togglePositions"``).

Example: Custom Tooltip Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: javascript

   const myTooltip = d3.select("body")
     .append("div")
     .attr("class", "tooltip")
     .style("opacity", 0);

   function customTooltip(event, d, html_content, isMouseOut) {
     if (isMouseOut) {
       myTooltip.transition().duration(200).style("opacity", 0);
     } else {
       myTooltip.html(html_content)
         .style("left", (event.pageX + 10) + "px")
         .style("top", (event.pageY - 28) + "px")
         .transition().duration(200).style("opacity", .9);
     }
   }

   createTreeMap("treeData.json", "#myTreeContainer", {
     tooltip: customTooltip
   }); 