# tree_plot.js Documentation

## Overview

`tree_plot.js` is a JavaScript library for visualizing transmission tree structures using D3.js. It supports interactive features like tooltips to see host attributes and color customization.

---

## Installation

1. **Include D3.js** in your HTML:
   ```html
   <script src="https://d3js.org/d3.v6.min.js"></script>
   ```
2. **Include `tree_plot.js`** in your HTML:
   ```html
   <script src="https://www.maths.usyd.edu.au/u/oscarf/tree_layout/tree_plot.js"></script>
   ```

---

## Usage

### 1. Prepare your HTML
Add a container `<div id="myTreeContainer"></div>` where the tree will be rendered:
```html
<div id="myTreeContainer"></div>
```

### 2. Prepare your data
Your data should be in JSON format, representing a tree structure. The library supports two formats:

#### **Format 1: Direct tree object**
```json
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
```

#### **Format 2: Tree object inside a 'tree' attribute**
```json
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
```
The library will automatically detect and use the correct format.

### 3. Call `createTreeMap`
```js
createTreeMap("treeData.json", "#myTreeContainer", {
  sampledColor: "#1f77b4", // Optional: color for sampled nodes
  unsampledColor: "#ff7f0e", // Optional: color for unsampled nodes
  tooltip: customTooltip, // Optional: custom tooltip function or d3 selection
  toggleSwitchId: "mySwitch" // Optional: custom toggle switch ID
});
```

---

## API Reference

### `createTreeMap(jsonFilename, containerDiv, options)`
- **`jsonFilename`** (string): Path to your tree data JSON file.
- **`containerDiv`** (string): CSS selector for the container div (e.g., `#myTreeContainer`).
- **`options`** (object, optional):
  - **`sampledColor`** (string): Color for sampled nodes (default: `"dodgerblue"`).
  - **`unsampledColor`** (string): Color for unsampled nodes (default: `"darkorange"`).
  - **`tooltip`** (function or d3 selection): Custom tooltip handler.
    - If a function, it will be called as `tooltip(event, d, html_content, isMouseOut)`.
    - If a d3 selection, it will be used as the tooltip div.
  - **`toggleSwitchId`** (string): ID of the toggle switch element (default: `"togglePositions"`).

#### Example: Custom Tooltip Function
```js
const myTooltip = d3.select("body")
  .append("div")
  .attr("id", "myCustomTooltip")
  .style("position", "absolute")
  .style("background", "#fff")
  .style("border", "1px solid #ccc")
  .style("padding", "8px")
  .style("border-radius", "4px")
  .style("pointer-events", "none")
  .style("opacity", 0);

function customTooltip(event, d, html_content, isMouseOut) {
  if (isMouseOut) {
    myTooltip.transition().duration(200).style("opacity", 0);
  } else {
    myTooltip.html(html_content)
      .style("left", (event.pageX + 15) + "px")
      .style("top", (event.pageY - 20) + "px")
      .transition().duration(200).style("opacity", 1);
  }
}
```

---

## Features
- **Interactive tooltips** on node hover
- **Customizable node colors** for sampled/unsampled nodes
- **Responsive resizing** on window resize
- **Dropdown for coloring by attributes** (if present in data)
- **Toggle for infection time layout** (if toggle switch is present in HTML)

---

## Related Libraries

### Transmission Models Library
For generating and sampling epidemiological transmission trees networks, you can use the [Transmission Models Library](https://github.com/oscarcapote/transmission_models) by Oscar Fajardo-Fontiveros. This Python library provides:
- MCMC sampling of transmission trees
- Genetic and location-based priors
- Epidemic network inference
- Support for different data types (genetic distance, location, timing)

This library is particularly useful for creating the tree data that can then be visualized using `tree_plot.js`.

---

## Layout Switch: Changing Between Classic and Infection Time Layouts

The library supports switching between two layouts: the classic layout and the infection time layout where hosts are placed according to their infection time. This is controlled by a toggle switch (checkbox) in your HTML.

### How to Set Up the Switch
Add this to your HTML:
```html
<label for="togglePositions">Use Infection Time Layout: </label>
<input type="checkbox" id="togglePositions">
```
- The switch must have the id `togglePositions` **if you use the default**, or the custom ID you specify in the `toggleSwitchId` option. The value `togglePositions` is just the default and can be changed as needed.
- The event listener is automatically set up by the library code.

### Customizing the Switch ID
You can use a custom ID for the toggle switch by passing it in the options:
```js
createTreeMap("treeData.json", "#myTreeContainer", { toggleSwitchId: "myCustomSwitch" });
```
Then use that ID in your HTML:
```html
<input type="checkbox" id="myCustomSwitch">
```

### How It Works
- When the switch is **unchecked**, the tree uses the classic layout (nodes positioned by their default y value).
- When the switch is **checked**, the tree uses the infection time layout (nodes positioned by their infection time, using the `infTimeScale`).
- The transition is animated for a smooth effect.

### What Happens Internally
- The function `positionNodesByInfectionTime(useInfectionTime, node, link, infTimeScale)` is called whenever the switch is toggled.
- The event listener is set up in the code:
```js
document.getElementById(toggleSwitchId).addEventListener('change', function() {
  positionNodesByInfectionTime(this.checked, node, link, infTimeScale);
});
```
- The plot updates automatically when the switch is toggled.

---

## Example
```html
<!-- HTML -->
<div id="myTreeContainer"></div>
<label for="togglePositions">Use Infection Time Layout: </label>
<input type="checkbox" id="togglePositions">

<script src="https://d3js.org/d3.v6.min.js"></script>
<script src="https://www.maths.usyd.edu.au/u/oscarf/tree_layout/tree_plot.js"></script>
<script>
  createTreeMap("treeData.json", "#myTreeContainer");
</script>
``` 