# Quick Start Guide - Tree Visualization Library

Get started with the Tree Visualization Library in just a few minutes.

## Prerequisites

- A web browser with JavaScript enabled
- D3.js (automatically loaded via CDN in examples)
- JSON data from the `transmission_models` package

## Step 1: Basic Setup

Create a simple HTML file with the following content:

```html
<!DOCTYPE html>
<html>
<head>
    <title>Tree Visualization</title>
    <script src="https://d3js.org/d3.v5.min.js"></script>
    <script src="https://www.maths.usyd.edu.au/u/oscarf/tree_layout/tree_plot.js"></script>
</head>
<body>
    <div id="tree-container"></div>
    
    <script>
        // Load and display your tree data
        fetch('path/to/your/model.json')
            .then(response => response.json())
            .then(data => {
                plot_Tree(data, "#tree-container");
            });
    </script>
</body>
</html>
```

## Step 2: Add Interactive Features

Enhance your visualization with a toggle switch for time-based layout:

```html
<!DOCTYPE html>
<html>
<head>
    <title>Interactive Tree Visualization</title>
    <script src="https://d3js.org/d3.v5.min.js"></script>
    <script src="https://www.maths.usyd.edu.au/u/oscarf/tree_layout/tree_plot.js"></script>
</head>
<body>
    <div>
        <label>Use Infection Time Layout: 
            <input type="checkbox" id="toggleSwitch">
        </label>
    </div>
    <div id="tree-container"></div>
    
    <script>
        fetch('path/to/your/model.json')
            .then(response => response.json())
            .then(data => {
                plot_Tree(data, "#tree-container");
            });
    </script>
</body>
</html>
```

## Step 3: Understanding Your Data

The library works with JSON files that have this structure:

```json
{
    "log_likelihood": 427.14,
    "tree": {
        "name": "Virtual_host",
        "index": -1,
        "Infection time": -34.34,
        "Sampled": false,
        "root": true,
        "children": [
            {
                "name": "hCoV-19/Australia/NSW1649/2021",
                "index": 0,
                "Infection time": -7.29,
                "Sampled": true,
                "Sampling time": 0,
                "root": false,
                "children": []
            }
        ]
    },
    "parameters": {
        "sampling_params": {...},
        "offspring_params": {...},
        "infection_params": {...}
    }
}
```

## Step 4: Key Features to Try

### 1. Hover Effects
- Move your mouse over any node to see it and its relatives highlight
- This helps identify transmission chains

### 2. Time Layout Toggle
- Check the toggle box to switch to infection time-based layout
- This shows the temporal progression of infections

### 3. Color Coding
- **Blue nodes**: Sampled hosts (with sequence data)
- **Orange nodes**: Unsampled hosts (inferred from transmission model)

### 4. Responsive Design
- Resize your browser window to see the tree adapt automatically

## Common Use Cases

### Display a Single Tree
```javascript
plot_Tree(data, "#container");
```

### Display Multiple Trees
```javascript
// Load multiple models
Promise.all([
    fetch('model1.json').then(r => r.json()),
    fetch('model2.json').then(r => r.json())
]).then(([data1, data2]) => {
    plot_Tree(data1, "#tree1");
    plot_Tree(data2, "#tree2");
});
```

### Custom Styling
```css
/* Customize node appearance */
.node circle {
    stroke: #333;
    stroke-width: 2px;
}

/* Customize link appearance */
.link {
    stroke: #666;
    stroke-width: 1.5px;
}
```

## Troubleshooting

### "d3 is not defined" Error
- Make sure D3.js is loaded before the tree visualization library
- Check that the D3.js CDN link is working

### "Container not found" Error
- Verify the CSS selector matches an existing element
- Ensure the element exists before calling `plot_Tree()`

### Tree Not Displaying
- Check browser console for JavaScript errors
- Verify your JSON data follows the expected format
- Ensure all nodes have required properties like "Infection time"

## Next Steps

- Read the [full documentation](tree_visualization.md) for advanced features
- Try the [interactive example](tree_visualization_example.html)
- Explore customization options for your specific needs

## Support

If you encounter issues:
1. Check the browser console for error messages
2. Verify your data format matches the expected structure
3. Ensure all dependencies are properly loaded
4. Test with the provided example files first 