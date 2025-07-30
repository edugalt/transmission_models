# Tree Visualization Library Documentation

This directory contains comprehensive documentation for the Tree Visualization Library, a JavaScript library for visualizing transmission tree data from the `transmission_models` Python package.

## Documentation Structure

### 📖 [index.md](index.md)
**Main documentation hub** - Overview, navigation, and quick reference for the entire library.

### 🚀 [quick_start.md](quick_start.md)
**Quick Start Guide** - Get up and running in minutes with basic examples and common use cases.

### 📚 [tree_visualization.md](tree_visualization.md)
**Complete Documentation** - Full API reference, detailed examples, customization options, and troubleshooting guide.

### 🎯 [tree_visualization_example.html](tree_visualization_example.html)
**Interactive Example** - A complete, working HTML page that demonstrates all library features with sample data.

## Getting Started

1. **New Users**: Start with [quick_start.md](quick_start.md) for a 5-minute introduction
2. **Interactive Demo**: Open [tree_visualization_example.html](tree_visualization_example.html) in your browser
3. **Full Reference**: Read [tree_visualization.md](tree_visualization.md) for complete documentation
4. **Navigation**: Use [index.md](index.md) as your documentation hub

## Library Overview

The Tree Visualization Library provides:

- **Interactive tree visualizations** using D3.js
- **Dual layout modes**: Hierarchical and infection time-based
- **Sampling status visualization**: Blue for sampled, orange for unsampled hosts
- **Interactive features**: Hover effects, responsive design, smooth transitions
- **Easy integration**: Simple API that works with JSON from transmission_models

## Quick Example

```html
<!DOCTYPE html>
<html>
<head>
    <script src="https://d3js.org/d3.v5.min.js"></script>
    <script src="https://www.maths.usyd.edu.au/u/oscarf/tree_layout/tree_plot.js"></script>
</head>
<body>
    <div id="tree-container"></div>
    
    <script>
        fetch('your_model.json')
            .then(response => response.json())
            .then(data => {
                plot_Tree(data, "#tree-container");
            });
    </script>
</body>
</html>
```

## Data Format

The library expects JSON data with this structure:

```json
{
    "log_likelihood": 427.14,
    "tree": {
        "name": "Virtual_host",
        "index": -1,
        "Infection time": -34.34,
        "Sampled": false,
        "root": true,
        "children": [...]
    },
    "parameters": {...}
}
```

## Features

### Core Features
- ✅ Interactive tree visualization
- ✅ Dual layout modes (hierarchical/time-based)
- ✅ Sampling status visualization
- ✅ Hover effects and highlighting
- ✅ Responsive design
- ✅ Smooth transitions

### Supported Model Types
- ✅ Simple transmission models
- ✅ Genetic models with distance data
- ✅ Location models with geographic data

### Browser Support
- ✅ Chrome 60+
- ✅ Firefox 55+
- ✅ Safari 12+
- ✅ Edge 79+

## Documentation Style

This documentation follows the style of the [BiMMSBM documentation](https://oscarcapote.github.io/BiMMSBM/index.html), providing:

- Clear navigation and structure
- Comprehensive examples
- Interactive demonstrations
- Troubleshooting guides
- API reference

## Contributing to Documentation

To improve the documentation:

1. **Report Issues**: If you find errors or unclear sections
2. **Suggest Improvements**: Propose better examples or explanations
3. **Add Examples**: Share your use cases and code samples
4. **Update Content**: Keep documentation current with library changes

## Library Information

- **Library URL**: https://www.maths.usyd.edu.au/u/oscarf/tree_layout/tree_plot.js
- **Dependencies**: D3.js v5+
- **License**: Provided as-is for use with transmission_models
- **Version**: 1.0.0

## Support

For questions and issues:

1. Check the [troubleshooting section](tree_visualization.md#troubleshooting)
2. Review the [data format requirements](tree_visualization.md#data-format-compatibility)
3. Test with the [interactive example](tree_visualization_example.html)
4. Ensure all dependencies are properly loaded

---

**Last Updated**: December 2024 