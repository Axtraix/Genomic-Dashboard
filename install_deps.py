import subprocess
import sys


required_libs = ["requests", "pandas", "pydantic", "pyarrow", "plotly", "numpy"]

print("Starting dependency check...")
for lib in required_libs:
    try:
        __import__(lib)
        print(f" -> {lib} is already installed.")
    except ImportError:
        print(f" -> Installing {lib}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", lib])

print("\nEnvironment setup complete!")