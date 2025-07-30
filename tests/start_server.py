#!/usr/bin/env python3
"""
Simple script to start a local web server for running the tree visualization tests.
This solves the CORS issue when trying to load local JSON files in the browser.
"""

import http.server
import socketserver
import os
import sys
import webbrowser
from pathlib import Path

def main():
    print("🌳 Tree Visualization Test Server")
    print("=" * 40)
    
    # Get the project root directory (parent of tests directory)
    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    
    print(f"Script location: {current_dir}")
    print(f"Project root: {project_root}")
    print(f"Current working directory: {os.getcwd()}")
    print()
    
    # Check if we're already in the project root
    if os.getcwd() != str(project_root):
        # Change to project root directory
        os.chdir(project_root)
        print(f"📁 Changed working directory to: {project_root}")
        print()
    
    # Verify that the example files exist
    example_files = [
        "examples/model_simple.json",
        "examples/model_genetic.json", 
        "examples/model_location.json"
    ]
    
    missing_files = []
    for file_path in example_files:
        if not Path(file_path).exists():
            missing_files.append(file_path)
    
    if missing_files:
        print("❌ Error: Missing example files:")
        for file_path in missing_files:
            print(f"   - {file_path}")
        print(f"\nCurrent directory: {os.getcwd()}")
        print("Expected project structure:")
        print("   transmission_models/")
        print("   ├── examples/")
        print("   │   ├── model_simple.json")
        print("   │   ├── model_genetic.json")
        print("   │   └── model_location.json")
        print("   ├── tests/")
        print("   │   └── start_server.py")
        print("   └── src/")
        print("\nPlease ensure you're running the script from the correct location.")
        print("The script should be run from the tests directory or the project root.")
        return
    
    print("✅ Example files verified:")
    for file_path in example_files:
        print(f"   • {file_path}")
    print()
    
    # Set up server with automatic port selection
    Handler = http.server.SimpleHTTPRequestHandler
    PORT = None
    
    # Try different ports if 8000 is in use
    for port in range(8000, 8010):
        try:
            with socketserver.TCPServer(("", port), Handler) as httpd:
                PORT = port
                break
        except OSError as e:
            if e.errno == 48:  # Address already in use
                continue
            else:
                raise e
    else:
        print("❌ Error: Could not find an available port between 8000-8009")
        return
    
    try:
        with socketserver.TCPServer(("", PORT), Handler) as httpd:
            print("🚀 Starting server...")
            print(f"Server started at: http://localhost:{PORT}")
            print()
            print("Available test files:")
            print(f"  • Test Suite: http://localhost:{PORT}/tests/test_suite.html")
            print(f"  • Simple Model: http://localhost:{PORT}/tests/test_simple_model.html")
            print(f"  • Genetic Model: http://localhost:{PORT}/tests/test_genetic_model.html")
            print(f"  • Location Model: http://localhost:{PORT}/tests/test_location_model.html")
            print()
            print("Press Ctrl+C to stop the server")
            print("=" * 60)
            
            # Try to open the test suite in the default browser
            try:
                webbrowser.open(f"http://localhost:{PORT}/tests/test_suite.html")
                print("✅ Opened test suite in your default browser")
            except:
                print("⚠️  Could not open browser automatically. Please open manually:")
                print(f"   http://localhost:{PORT}/tests/test_suite.html")
            
            print()
            httpd.serve_forever()
            
    except KeyboardInterrupt:
        print("\n🛑 Server stopped by user")
    except OSError as e:
        if e.errno == 48:  # Address already in use
            print(f"❌ Port {PORT} is already in use. Try a different port:")
            print(f"   python3 -m http.server 8001")
        else:
            print(f"❌ Error starting server: {e}")
    except Exception as e:
        print(f"❌ Unexpected error: {e}")

if __name__ == "__main__":
    main() 