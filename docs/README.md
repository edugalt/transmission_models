# Transmission Models Documentation

This directory contains the Sphinx documentation for the transmission_models package.

## Building the Documentation

### Prerequisites

Install the required dependencies:

```bash
pip install -r requirements.txt
```

### Build Commands

1. **Build HTML documentation:**
   ```bash
   make html
   ```

2. **Build PDF documentation:**
   ```bash
   make latexpdf
   ```

3. **Clean build directory:**
   ```bash
   make clean
   ```

4. **View available build targets:**
   ```bash
   make help
   ```

### Viewing the Documentation

After building with `make html`, you can view the documentation by opening:
```
_build/html/index.html
```

## Documentation Structure

- `conf.py`: Sphinx configuration file
- `index.rst`: Main documentation page
- `api.rst`: API reference documentation
- `requirements.txt`: Python dependencies for documentation
- `Makefile`: Build automation

## Configuration

The documentation is configured to use:
- NumPy-style docstrings
- Read the Docs theme
- Automatic API documentation generation
- Intersphinx links to Python, NumPy, SciPy, NetworkX, and Matplotlib

## Adding New Documentation

1. Add new RST files to the documentation
2. Update `index.rst` to include new pages in the toctree
3. Update `api.rst` to include new modules
4. Rebuild the documentation with `make html` 