# Tree Visualization Library Tests

This directory contains comprehensive HTML test files for the Tree Visualization Library, designed to validate functionality with different types of transmission model data.

## Test Files

### 🧪 Individual Model Tests

#### [test_simple_model.html](test_simple_model.html)
**Simple Model Test** - Tests basic tree visualization with a simple transmission model.
- Tests data loading and structure validation
- Validates tree rendering and node count
- Checks sampling status visualization
- Tests interactive features (toggle switch, hover effects)

#### [test_genetic_model.html](test_genetic_model.html)
**Genetic Model Test** - Tests tree visualization with genetic distance data.
- All simple model tests plus genetic-specific tests
- Validates genetic prior calculations
- Tests visualization with genetic model data
- Ensures proper handling of genetic distance information

#### [test_location_model.html](test_location_model.html)
**Location Model Test** - Tests tree visualization with geographic location data.
- All simple model tests plus location-specific tests
- Validates spatial constraints and location data
- Tests visualization with location model data
- Ensures proper handling of geographic information

### 🎯 [test_suite.html](test_suite.html)
**Comprehensive Test Suite** - Runs all model tests together with summary reporting.
- Tests all three model types simultaneously
- Provides overall test statistics and summary
- Interactive test controls and result export
- Side-by-side comparison of different model visualizations

## Running the Tests

### Prerequisites
- Modern web browser (Chrome, Firefox, Safari, Edge)
- Internet connection (for D3.js CDN)
- Access to the example JSON files in `../examples/`

### ⚠️ Important: CORS Issue Solution

**Problem**: If you see "Cross-Origin Request Blocked" errors, this is because browsers block loading local files directly for security reasons.

**Solution**: Use a local web server instead of opening files directly.

### Quick Start (Recommended)

#### Option 1: Use the provided Python script
```bash
# Run from anywhere in the project
python3 tests/start_server.py

# Or run from the tests directory
cd tests
python3 start_server.py
```

This will:
- ✅ Automatically find the project root directory
- ✅ Start a local web server on port 8000
- ✅ Verify all example files exist
- ✅ Automatically open the test suite in your browser
- ✅ Show you all available test URLs

#### Option 2: Manual server setup
```bash
# Navigate to the project root
cd /home/oscarf/Documents/transmission_models

# Start a Python HTTP server
python3 -m http.server 8000

# Then open in your browser:
# http://localhost:8000/tests/test_suite.html
```

#### Option 3: Alternative servers
```bash
# Node.js (if you have Node.js installed)
npx http-server -p 8000

# PHP (if you have PHP installed)
php -S localhost:8000
```

### Test URLs (after starting server)
- **Test Suite**: http://localhost:8000/tests/test_suite.html
- **Simple Model**: http://localhost:8000/tests/test_simple_model.html
- **Genetic Model**: http://localhost:8000/tests/test_genetic_model.html
- **Location Model**: http://localhost:8000/tests/test_location_model.html

### Individual Model Testing
```bash
# Open specific model tests
open tests/test_simple_model.html
open tests/test_genetic_model.html
open tests/test_location_model.html
```

### Comprehensive Testing
```bash
# Open the full test suite
open tests/test_suite.html
```

## Test Features

### ✅ Data Validation Tests
- **Data Loading**: Verifies JSON files can be loaded successfully
- **Data Structure**: Validates required fields are present
- **Model-Specific Fields**: Checks for genetic prior, location data, etc.
- **Tree Structure**: Ensures proper node hierarchy and relationships

### ✅ Visualization Tests
- **Rendering**: Confirms SVG elements are created
- **Node Count**: Validates correct number of nodes displayed
- **Interactive Features**: Tests hover effects and toggle functionality
- **Responsive Design**: Ensures visualization adapts to container size

### ✅ Model-Specific Tests
- **Simple Model**: Basic transmission tree functionality
- **Genetic Model**: Genetic prior calculations and distance data
- **Location Model**: Geographic constraints and spatial data

### ✅ User Interface Tests
- **Toggle Switch**: Tests infection time layout switching
- **Error Handling**: Validates proper error messages
- **Loading States**: Ensures appropriate loading indicators
- **Result Display**: Tests test result formatting and display

## Test Results

### Pass/Fail Indicators
- 🟢 **Pass**: Test completed successfully
- 🔴 **Fail**: Test failed with error details
- 🟡 **Warning**: Test completed with warnings

### Summary Statistics
- **Total Tests**: Number of tests run
- **Passed**: Number of successful tests
- **Failed**: Number of failed tests
- **Warnings**: Number of tests with warnings

### Export Functionality
- **JSON Export**: Download test results as JSON file
- **Timestamp**: Results include test execution timestamp
- **Detailed Results**: Full test output with model-specific information

## Test Data

The tests use the following example files from the `transmission_models` package:

### Model Files
- `../examples/model_simple.json` - Basic transmission model
- `../examples/model_genetic.json` - Model with genetic distance data
- `../examples/model_location.json` - Model with geographic location data

### Data Structure Validation
Each test validates the expected JSON structure:
```json
{
    "log_likelihood": number,
    "tree": {
        "name": string,
        "index": number,
        "Infection time": number,
        "Sampled": boolean,
        "root": boolean,
        "children": [...]
    },
    "parameters": {...},
    "genetic_prior": number,  // Optional (genetic model)
    "location_prior": number  // Optional (location model)
}
```

## Troubleshooting

### Common Issues

#### "Cross-Origin Request Blocked" (CORS Error)
- **Problem**: Browser blocks loading local files directly
- **Solution**: Use a local web server (see Quick Start section above)
- **Quick Fix**: Run `python3 start_server.py` in the tests directory
- **Alternative**: Use `python3 -m http.server 8000` in project root

#### "Failed to load model data"
- **Solution**: Ensure example JSON files exist in `../examples/`
- **Check**: File paths and permissions
- **Verify**: JSON files contain valid data
- **Note**: Must access via `http://localhost:8000` not `file://`

#### "Visualization not rendered"
- **Solution**: Check browser console for JavaScript errors
- **Verify**: D3.js is loaded successfully
- **Check**: Tree visualization library is accessible
- **Common Cause**: CORS error preventing data loading

#### "Toggle switch not working"
- **Solution**: Ensure toggle element has correct ID (`toggleSwitch`)
- **Check**: Event listeners are properly attached
- **Verify**: Library supports toggle functionality

#### "Port 8000 is already in use"
- **Solution**: Use a different port
- **Try**: `python3 -m http.server 8001`
- **Then**: Open `http://localhost:8001/tests/test_suite.html`

### Debug Mode
Enable browser developer tools to see:
- Console errors and warnings
- Network requests for JSON files
- DOM structure of visualizations
- JavaScript execution flow

## Continuous Integration

### Automated Testing
These HTML tests can be integrated into CI/CD pipelines using:
- **Selenium WebDriver** for automated browser testing
- **Puppeteer** for headless browser automation
- **Playwright** for cross-browser testing

### Test Reporting
- **JSON Export**: Machine-readable test results
- **Console Output**: Browser console logging
- **Visual Feedback**: Real-time test status updates

## Contributing to Tests

### Adding New Tests
1. **Create test file** following existing naming convention
2. **Include test functions** for specific functionality
3. **Add to test suite** if applicable
4. **Update documentation** with new test details

### Test Best Practices
- **Isolated Tests**: Each test should be independent
- **Clear Messages**: Provide descriptive error messages
- **Comprehensive Coverage**: Test both success and failure cases
- **Performance**: Ensure tests run efficiently

### Test Maintenance
- **Regular Updates**: Keep tests current with library changes
- **Data Validation**: Ensure test data remains valid
- **Browser Compatibility**: Test across different browsers
- **Performance Monitoring**: Track test execution times

## Support

For test-related issues:
1. **Check browser console** for error messages
2. **Verify file paths** and data accessibility
3. **Test with different browsers** to isolate issues
4. **Review test documentation** for specific guidance

---

**Last Updated**: December 2024 